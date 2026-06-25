# Conformer clustering pipeline -------------------------------------------------------

# Reference handling
# -----

# Normalize apo/holo reference stem: blank/whitespace → nothing, else unchanged.
function _normalize_reference(ref::String)::Union{Nothing,String}
    stripped = strip(ref)
    return isempty(stripped) ? nothing : stripped
end

# Detect reference mode from normalized apo/holo stems.
function _reference_mode(apo::Union{Nothing,String}, holo::Union{Nothing,String})::Int
    if holo !== nothing && apo === nothing
        error("holo reference given without apo reference; use 0-ref or 1-ref mode.")
    end
    count = (apo !== nothing) + (holo !== nothing)
    return count
end

# Conformer discovery
# -----

# Directory name that holds the AlphaConformer prediction structures in the package's output
# tree (`.../predictions/sequences/models/`). Only files under such a directory are treated as
# predictions; this keeps input templates (`templates_complete/`, `templates_adaptative/`) and
# search results (`fullpdb_results/`, `target_db_results/`) out of the ensemble.
const _MODELS_DIR = "models"

# Discover the AlphaConformer prediction structures of one system below `data_root/<system>`.
#
# Recurses the system directory and keeps only structure files whose immediate parent directory
# is a `models/` folder (the AlphaConformers prediction location). Files under template or
# search-result folders are skipped, so templates are never mistaken for predictions. Returns a
# vector of `(name, path)` tuples sorted by `name`, where `name` is the file path relative to
# the system directory (stable across machines) and `path` is the absolute path. Raises a clear
# error if the system directory is missing or holds no prediction structures.
function _discover_conformers(data_root::AbstractString, system::AbstractString)
    system_dir = joinpath(data_root, system)
    isdir(system_dir) ||
        error("no conformer directory found for system '$system' at '$system_dir'")

    found = Tuple{String,String}[]
    for (root, _dirs, files) in walkdir(system_dir)
        basename(root) == _MODELS_DIR || continue
        for f in files
            ext = lowercase(splitext(f)[2])
            if ext == ".pdb" || ext == ".cif"
                path = joinpath(root, f)
                name = relpath(path, system_dir)
                push!(found, (name, abspath(path)))
            end
        end
    end

    isempty(found) && error(
        "no .pdb/.cif prediction structures found for system '$system' under '$system_dir'; " *
        "predictions are read only from `$_MODELS_DIR/` folders (the AlphaConformers " *
        "prediction location), so input templates are not treated as conformers.",
    )
    sort!(found; by = first)
    return found
end

# CA parsing
# -----

# Parse the C-alpha trace of one structure file.
#
# Returns `(resids, coords)` where `resids` is the vector of residue-number strings that carry
# a CA atom and `coords` is the matching `length(resids)×3` matrix of CA coordinates. Reads the
# whole file, following the single-chain-per-file convention.
function _parse_ca(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    residues = if ext == ".cif"
        MIToS.PDB.read_file(path, MIToS.PDB.MMCIFFile)
    else
        MIToS.PDB.read_file(path, MIToS.PDB.PDBFile)
    end

    resids = String[]
    rows = Vector{NTuple{3,Float64}}()
    for res in residues
        cas = MIToS.PDB.findatoms(res, "CA")
        isempty(cas) && continue
        c = res.atoms[first(cas)].coordinates
        push!(resids, res.id.number)
        push!(rows, (c.x, c.y, c.z))
    end

    coords = Matrix{Float64}(undef, length(rows), 3)
    for (i, r) in enumerate(rows)
        coords[i, 1], coords[i, 2], coords[i, 3] = r
    end
    return resids, coords
end

# Common CA positions shared by every conformer, keeping the order of the first conformer.
#
# `resids_list` is a vector of residue-number-string vectors, one per conformer. Returns the
# ordered intersection. In reference-free mode every prediction shares one sequence, so this
# is normally the full set; conformers with a few missing residues still line up cleanly.
function _common_ca_positions(resids_list::AbstractVector)
    isempty(resids_list) && error("no conformers to find common CA positions for")
    common = Set(resids_list[1])
    for resids in Iterators.drop(resids_list, 1)
        intersect!(common, Set(resids))
    end
    isempty(common) &&
        error("conformers share no common CA residue positions; cannot cluster")
    return [r for r in resids_list[1] if r in common]
end

# Restrict one conformer's CA coordinates to the common positions, in `common` order.
function _restrict_to_common(
    resids::AbstractVector,
    coords::AbstractMatrix,
    common::AbstractVector,
)
    index = Dict(r => i for (i, r) in enumerate(resids))
    out = Matrix{Float64}(undef, length(common), 3)
    for (i, r) in enumerate(common)
        row = index[r]
        out[i, 1], out[i, 2], out[i, 3] = coords[row, 1], coords[row, 2], coords[row, 3]
    end
    return out
end

# Superposition (folded from scripts/conformer_clustering/superimpose_helper.jl)
# -----

# Build a vector of CA-only `PDBResidue`s from an `L×3` coordinate matrix.
function _make_ca_residues(coords::AbstractMatrix)
    n = size(coords, 1)
    residues = Vector{MIToS.PDB.PDBResidue}(undef, n)
    for i = 1:n
        atom = MIToS.PDB.PDBAtom(
            MIToS.PDB.Coordinates(coords[i, 1], coords[i, 2], coords[i, 3]),
            "CA",
            "C",
            1.0,
            "0.00",
            " ",
            " ",
        )
        res_id = MIToS.PDB.PDBResidueIdentifier(" ", string(i), "ALA", "ATOM", "1", "A")
        residues[i] = MIToS.PDB.PDBResidue(res_id, [atom])
    end
    return residues
end

# Extract the `L×3` coordinate matrix back from a vector of CA `PDBResidue`s.
function _residues_to_matrix(residues::AbstractVector)
    n = length(residues)
    out = Matrix{Float64}(undef, n, 3)
    for (i, r) in enumerate(residues)
        c = r.atoms[1].coordinates
        out[i, 1], out[i, 2], out[i, 3] = c.x, c.y, c.z
    end
    return out
end

# Superimpose every conformer onto the first one (all share the same `L` CA positions).
#
# `coords_list` is a vector of `L×3` matrices. Returns a matching vector of `L×3` matrices in
# the first conformer's centered frame. Uses the package's MIToS superposition as the single
# source of truth for structural alignment.
function _superimpose_to_first(coords_list::AbstractVector)
    ref = _make_ca_residues(coords_list[1])
    ref_centered, _, _ = MIToS.PDB.superimpose(ref, ref)
    aligned = Vector{Matrix{Float64}}(undef, length(coords_list))
    aligned[1] = _residues_to_matrix(ref_centered)
    for i = 2:length(coords_list)
        _, bsuper, _ = MIToS.PDB.superimpose(ref, _make_ca_residues(coords_list[i]))
        aligned[i] = _residues_to_matrix(bsuper)
    end
    return aligned
end

# Build the `(3L)×N` feature matrix (one observation per column) from aligned coordinates.
function _feature_matrix(aligned::AbstractVector)
    n = length(aligned)
    d = length(aligned[1])
    features = Matrix{Float64}(undef, d, n)
    for i = 1:n
        features[:, i] = vec(aligned[i])
    end
    return features
end

# Reference RMSD
# -----

# Resolve an apo/holo reference stem (`<PDBID>_<CHAIN>`) to a file under `refs_dir`.
#
# Tries, across `.pdb` then `.cif`: the exact stem first, then a case-insensitive
# lowercase-id fallback that lowercases the PDB id but keeps the chain letter
# (`1AKZ_A` → `1akz_A`). Returns the first existing path. Raises a clear error when the stem
# resolves to no file. No chain argument: each reference file is read whole.
function _resolve_reference(refs_dir::AbstractString, stem::AbstractString)
    candidates = String[String(stem)]
    if occursin('_', stem)
        pdb_id, chain_letter = rsplit(String(stem), '_'; limit = 2)
        push!(candidates, "$(lowercase(pdb_id))_$(chain_letter)")
    end
    for cand in candidates, ext in (".pdb", ".cif")
        path = joinpath(refs_dir, cand * ext)
        isfile(path) && return path
    end
    tried = join([cand * ext for cand in candidates for ext in (".pdb", ".cif")], ", ")
    error("could not resolve reference stem '$stem' under '$refs_dir'; tried: $tried")
end

# Map a few common modified residues back to their standard parents so sequence extraction
# and alignment treat them as the canonical residue. Mirrors the reference RMSD script.
const _MODIFIED_AA = Dict(
    "MSE" => "MET",  # selenomethionine
    "SEC" => "CYS",  # selenocysteine
    "SEP" => "SER",  # phosphoserine
    "TPO" => "THR",  # phosphothreonine
    "PTR" => "TYR",  # phosphotyrosine
    "HYP" => "PRO",  # hydroxyproline
    "MLY" => "LYS",  # N-methyllysine
    "CSO" => "CYS",  # S-hydroxycysteine
    "PCA" => "GLN",  # pyroglutamate
)

# Rename modified residues in place to their standard parent.
function _normalize_modified_residues!(residues)
    for r in residues
        std = get(_MODIFIED_AA, r.id.name, nothing)
        std === nothing && continue
        old = r.id
        r.id = MIToS.PDB.PDBResidueIdentifier(
            old.PDBe_number,
            old.number,
            std,
            "ATOM",
            old.model,
            old.chain,
        )
    end
    return residues
end

# Read one structure file into its CA-bearing `PDBResidue`s for sequence-aware RMSD.
#
# Reads the whole file (`.cif` → MMCIF, otherwise PDB), normalizes modified residues, and
# keeps only residues that carry a CA atom. No chain argument: the single-chain-per-file
# convention means the whole file is one chain.
function _read_ca_residues(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    residues = if ext == ".cif"
        MIToS.PDB.read_file(path, MIToS.PDB.MMCIFFile)
    else
        MIToS.PDB.read_file(path, MIToS.PDB.PDBFile)
    end
    _normalize_modified_residues!(residues)
    filter!(r -> !isempty(MIToS.PDB.findatoms(r, "CA")), residues)
    return residues
end

# Sequence-aware RMSD from one reference to every conformer.
#
# Reuses `structural_alignment` (the package's single source of truth for structural
# distance) per conformer; an alignment that fails (`nothing`) yields `NaN` so the
# conformer drops out of the reference's cluster mean.
function _reference_rmsds(ref_residues::AbstractVector, conformer_residues::AbstractVector)
    rmsds = Vector{Float64}(undef, length(conformer_residues))
    for (i, conf) in enumerate(conformer_residues)
        result = structural_alignment(ref_residues, conf)
        rmsds[i] = result === nothing ? NaN : result[4]
    end
    return rmsds
end

# Assign a reference to the cluster whose conformers have the lowest mean RMSD to it.
#
# `rmsds` is the per-conformer RMSD to the reference and `labels` the matching cluster
# labels. NaN RMSDs are ignored, and clusters whose RMSDs are all NaN are skipped. Returns
# `(cluster, mean_rmsd)`. Raises when no cluster has a finite RMSD.
function _assign_reference_cluster(rmsds::AbstractVector, labels::AbstractVector)
    best_cluster = -1
    best_mean = Inf
    for lbl in sort(unique(labels))
        idx = findall(==(lbl), labels)
        vals = filter(!isnan, rmsds[idx])
        isempty(vals) && continue
        m = Statistics.mean(vals)
        if m < best_mean
            best_mean = m
            best_cluster = lbl
        end
    end
    best_cluster == -1 &&
        error("could not assign reference to any cluster; all RMSDs are NaN")
    return best_cluster, best_mean
end

# Write the apo/holo RMSD table and the reference cluster-assignment table.
#
# `refs` is a vector of `(role, stem, rmsds)` tuples (`role` is `"apo"` or `"holo"`). Writes
# `dist_external_rmsds.csv` (one row per conformer, a column per reference role) and
# `reference_clusters.csv` (one row per reference), both at the per-system output root.
# Returns both tables.
function _write_reference_tables(
    system_dir::AbstractString,
    names::AbstractVector,
    labels::AbstractVector,
    refs::AbstractVector,
)
    rmsd_table = DataFrames.DataFrame(:Name => collect(names))
    assignment_rows = NamedTuple[]
    for (role, stem, rmsds) in refs
        col = role == "apo" ? :RMSD_Apo : :RMSD_Holo
        rmsd_table[!, col] = rmsds
        cluster, mean_rmsd = _assign_reference_cluster(rmsds, labels)
        push!(
            assignment_rows,
            (
                Reference = role,
                Stem = stem,
                Cluster = cluster,
                Mean_RMSD_Angstrom = mean_rmsd,
            ),
        )
    end

    mkpath(system_dir)
    CSV.write(joinpath(system_dir, "dist_external_rmsds.csv"), rmsd_table)

    reference_clusters = DataFrames.DataFrame(assignment_rows)
    CSV.write(joinpath(system_dir, "reference_clusters.csv"), reference_clusters)

    return rmsd_table, reference_clusters
end

# KMeans clustering
# -----

# Run seeded KMeans on the `(3L)×N` feature matrix.
#
# Returns 0-based labels (`0:k-1`) to match the Python pipeline's output contract. The RNG is
# seeded from `seed` so the same inputs reproduce the same clustering.
function _kmeans_cluster(features::AbstractMatrix, k::Integer; seed::Integer)
    rng = Random.MersenneTwister(seed)
    result = Clustering.kmeans(features, k; rng = rng)
    return Int.(Clustering.assignments(result)) .- 1
end

# Mean intra-cluster RMSD to the cluster centroid, per cluster label.
#
# `features` is the `(3L)×N` aligned-coordinate matrix, `labels` the 0-based assignments and
# `n_atoms` the number of CA positions. Mirrors the Python lean-mode centroid metric; singleton
# clusters get `0.0`.
function _intra_cluster_rmsd(
    features::AbstractMatrix,
    labels::AbstractVector,
    n_atoms::Integer,
)
    rmsd = Dict{Int,Float64}()
    for lbl in unique(labels)
        idx = findall(==(lbl), labels)
        if length(idx) < 2
            rmsd[lbl] = 0.0
            continue
        end
        members = features[:, idx]
        centroid = vec(Statistics.mean(members; dims = 2))
        per_member =
            [sqrt(sum(abs2, members[:, j] .- centroid) / n_atoms) for j in axes(members, 2)]
        rmsd[lbl] = Statistics.mean(per_member)
    end
    return rmsd
end

# Output tables
# -----

# Write the clustering and intra-cluster RMSD tables under `kmeans_dir` and return them.
#
# `aligned_clustering_results.csv` has columns `Name` and `KMeans_K<k>` (labels `0:k-1`).
# `aligned_cluster_rmsd.csv` has columns `Method`, `Cluster`, `Size`, `Mean_RMSD_Angstrom`,
# one row per occupied cluster.
function _write_kmeans_tables(
    kmeans_dir::AbstractString,
    names::AbstractVector,
    labels::AbstractVector,
    k::Integer,
    rmsd::AbstractDict,
)
    mkpath(kmeans_dir)
    method = "KMeans_K$k"

    clustering = DataFrames.DataFrame(:Name => collect(names), Symbol(method) => labels)
    CSV.write(joinpath(kmeans_dir, "aligned_clustering_results.csv"), clustering)

    occupied = sort(collect(keys(rmsd)))
    cluster_rmsd = DataFrames.DataFrame(
        Method = fill(method, length(occupied)),
        Cluster = occupied,
        Size = [count(==(lbl), labels) for lbl in occupied],
        Mean_RMSD_Angstrom = [rmsd[lbl] for lbl in occupied],
    )
    CSV.write(joinpath(kmeans_dir, "aligned_cluster_rmsd.csv"), cluster_rmsd)

    return clustering, cluster_rmsd
end

# DeepAccNet score filter
# -----

# Align a DeepAccNet score table to the discovered conformer names.
#
# `score_table` is either a `DataFrame` or a path to a CSV read with `CSV.read`. It must hold a
# `Name` column and the `score_col` score column. Names are matched against `names`; conformers
# absent from the table are dropped (a Python-like inner join), so the filter operates only on
# scored conformers. Returns `(scored_idx, scores)`, where `scored_idx` are the positions in
# `names` that have a score (in `names` order) and `scores` the matching score values. Raises an
# error if a required column is missing or if no conformer matches at all.
function _align_scores(
    score_table,
    names::AbstractVector;
    score_col::AbstractString = "cb-lddt",
)
    table =
        score_table isa DataFrames.AbstractDataFrame ? score_table :
        CSV.read(score_table, DataFrames.DataFrame)

    hasproperty(table, :Name) || error("score table is missing a `Name` column")
    hasproperty(table, Symbol(score_col)) ||
        error("score table is missing the score column `$score_col`")

    by_name = Dict(string(n) => s for (n, s) in zip(table.Name, table[!, score_col]))

    scored_idx = Int[]
    scores = Float64[]
    for (i, nm) in enumerate(names)
        s = get(by_name, string(nm), nothing)
        s === nothing && continue
        push!(scored_idx, i)
        push!(scores, Float64(s))
    end

    isempty(scored_idx) &&
        error("score table matched none of the conformers by `Name`; check the score table")
    return scored_idx, scores
end

# Determine which clusters survive a single relative/fraction cut.
#
# `labels` and `scores` are the scored subset (one entry per scored conformer). The relative
# cutoff is the `rel_pct` percentile of all scores; a member is "removed" when its score is
# strictly below the cutoff. A cluster survives when at most `frac_pct`% of its members are
# removed. Returns the sorted vector of surviving labels. Deterministic given the inputs.
function _surviving_clusters(
    labels::AbstractVector,
    scores::AbstractVector,
    rel_pct::Real,
    frac_pct::Real,
)
    cutoff = Statistics.quantile(scores, rel_pct / 100)
    below = scores .< cutoff
    frac = frac_pct / 100

    surviving = Int[]
    for lbl in sort(unique(labels))
        in_cluster = labels .== lbl
        n_total = count(in_cluster)
        frac_removed = n_total == 0 ? 1.0 : count(below .& in_cluster) / n_total
        frac_removed <= frac && push!(surviving, lbl)
    end
    return surviving
end

# Write the surviving-cluster table for one relative/fraction cut under
# `deepaccnet_dir/surviving/`.
#
# `names`, `labels` and `scores` are the scored subset. Computes the surviving clusters for the
# `rel`/`frac` cut and writes the scored conformers in those clusters to
# `surviving_rel<rel>_frac<frac>.csv` with columns `Name` and `KMeans_K<k>` (matching the kmeans
# clustering table). Only this configured cut is produced. Returns the sorted vector of
# surviving cluster labels.
function _write_surviving_table(
    deepaccnet_dir::AbstractString,
    names::AbstractVector,
    labels::AbstractVector,
    scores::AbstractVector,
    k::Integer,
    rel::Real,
    frac::Real,
)
    surv_dir = joinpath(deepaccnet_dir, "surviving")
    mkpath(surv_dir)
    method = "KMeans_K$k"

    surviving = _surviving_clusters(labels, scores, rel, frac)
    keep = Set(surviving)
    idx = findall(lbl -> lbl in keep, labels)
    table = DataFrames.DataFrame(
        :Name => collect(names[idx]),
        Symbol(method) => collect(labels[idx]),
    )
    CSV.write(joinpath(surv_dir, "surviving_rel$(rel)_frac$(frac).csv"), table)
    return surviving
end

# Agglomerative sub-clustering
# -----

# Full `N×N` pairwise CA-RMSD matrix for one set of conformers.
#
# `coords_list` is a vector of `L×3` CA-coordinate matrices (the members of one KMeans
# cluster). Each pair is optimally superimposed with the package's MIToS superposition (the
# single source of truth for structural distance), whose RMSD is rotation/translation
# invariant, so passing the already-aligned coordinates is correct. Returns a symmetric matrix
# with a zero diagonal.
function _pairwise_rmsd_matrix(coords_list::AbstractVector)
    n = length(coords_list)
    dist = zeros(Float64, n, n)
    residues = [_make_ca_residues(c) for c in coords_list]
    for i = 1:n, j = (i+1):n
        _, _, rmsd = MIToS.PDB.superimpose(residues[i], residues[j])
        dist[i, j] = rmsd
        dist[j, i] = rmsd
    end
    return dist
end

# Cut an agglomerative tree of one cluster into sub-clusters at a distance threshold.
#
# `dist` is the `N×N` pairwise RMSD matrix, `threshold` the cut height in Angstrom and
# `linkage` the linkage method (e.g. `:average`). Returns 1-based sub-cluster labels, one per
# member, in the input order. A cluster of fewer than two members cannot form a tree, so it
# becomes its own single sub-cluster (label `1`); an empty cluster returns no labels. The cut
# is deterministic given `dist`, `threshold` and `linkage`.
function _agglomerative_subcluster(dist::AbstractMatrix, threshold::Real, linkage::Symbol)
    n = size(dist, 1)
    n == 0 && return Int[]
    n == 1 && return [1]
    tree = Clustering.hclust(dist; linkage = linkage)
    return Int.(Clustering.cutree(tree; h = threshold))
end

# Build the per-conformer agglomerative sub-cluster assignment table.
#
# Columns: `Name`, `KMeans_K<k>` (the 0-based KMeans label, matching the kmeans table),
# `Sub_Cluster` (the 1-based label within the KMeans cluster) and `Mini_Cluster` (the combined
# string `c<label>_s<sub>` that uniquely names each mini-cluster across the ensemble). One row
# per conformer, in `names` order. Returns the table without writing it.
function _hierarchical_table(
    names::AbstractVector,
    kmeans_labels::AbstractVector,
    sub_labels::AbstractVector,
    k::Integer,
)
    mini = ["c$(lbl)_s$(sub)" for (lbl, sub) in zip(kmeans_labels, sub_labels)]
    return DataFrames.DataFrame(
        :Name => collect(names),
        Symbol("KMeans_K$k") => collect(kmeans_labels),
        :Sub_Cluster => collect(sub_labels),
        :Mini_Cluster => mini,
    )
end

# Write one KMeans cluster's membership table into its `kmeans<C>/` folder.
#
# Writes `members.csv`, this cluster's members with their sub-cluster labels (same schema as the
# flat top-level table: `Name`, `KMeans_K<k>`, `Sub_Cluster`, `Mini_Cluster`). Returns the
# cluster folder path.
function _write_cluster_outputs(
    cluster_dir::AbstractString,
    names::AbstractVector,
    kmeans_labels::AbstractVector,
    sub_labels::AbstractVector,
    k::Integer,
)
    mkpath(cluster_dir)
    members = _hierarchical_table(names, kmeans_labels, sub_labels, k)
    CSV.write(joinpath(cluster_dir, "members.csv"), members)
    return cluster_dir
end

# Dimensionality reduction (stdlib linear algebra)
# -----

# Principal-component projection of a `D×N` feature matrix onto `ndims` axes.
#
# Columns are observations. The columns are centered, then projected onto the leading left
# singular vectors via an SVD. Returns an `ndims×N` matrix of scores (one column per
# observation). When the data spans fewer than `ndims` directions the missing rows are zero,
# so the result is always `ndims×N`. Deterministic given the input.
function _pca(features::AbstractMatrix, ndims::Integer = 2)
    n = size(features, 2)
    mu = vec(Statistics.mean(features; dims = 2))
    centered = features .- mu
    decomposition = LinearAlgebra.svd(centered)
    k = min(ndims, size(decomposition.U, 2))
    scores = decomposition.U[:, 1:k]' * centered
    k < ndims && (scores = vcat(scores, zeros(ndims - k, n)))
    return Matrix(scores)
end

# Classical multidimensional scaling (PCoA) of an `N×N` distance matrix into `ndims`.
#
# Double-centers the squared distances and eigendecomposes the resulting Gram matrix, placing
# each point on the axes of the leading positive eigenvalues. Returns an `ndims×N` coordinate
# matrix (one column per point); axes with a non-positive eigenvalue stay zero. Recovers the
# original geometry up to rotation/reflection. Deterministic given the input.
function _classical_mds(dist::AbstractMatrix, ndims::Integer = 2)
    n = size(dist, 1)
    n == 0 && return zeros(ndims, 0)
    squared = dist .^ 2
    centering = LinearAlgebra.I - fill(1 / n, n, n)
    gram = -0.5 .* (centering * squared * centering)
    gram = (gram .+ gram') ./ 2
    decomposition = LinearAlgebra.eigen(LinearAlgebra.Symmetric(gram))
    order = sortperm(decomposition.values; rev = true)
    coords = zeros(ndims, n)
    for d = 1:min(ndims, n)
        idx = order[d]
        value = decomposition.values[idx]
        value <= 0 && continue
        coords[d, :] = decomposition.vectors[:, idx] .* sqrt(value)
    end
    return coords
end

# Full `N×N` euclidean distance matrix between the columns of a `D×N` feature matrix.
#
# Used to feed classical MDS for the embedding views; symmetric with a zero diagonal.
function _feature_distance_matrix(features::AbstractMatrix)
    n = size(features, 2)
    dist = zeros(Float64, n, n)
    for i = 1:n, j = (i+1):n
        d = sqrt(sum(abs2, view(features, :, i) .- view(features, :, j)))
        dist[i, j] = d
        dist[j, i] = d
    end
    return dist
end

# Plotting (kept separate from the algorithmic functions)
# -----

# Every figure is saved at this device resolution and figure size, matching the legacy Python
# pipeline (600 dpi, 13×8 inches; Plots sizes in pixels at its 100 px/inch convention). The
# margin keeps the axis titles from being clipped at the figure edge (Plots has no automatic
# tight layout).
const _FIG_DPI = 600
const _FIG_SIZE = (1300, 800)
const _FIG_MARGIN = 10plt.mm

# Colors a `-1` noise label as light gray, mirroring the Python `_cluster_color_map`.
const _NOISE_COLOR = plt.RGBA(0.7, 0.7, 0.7, 1.0)

# Map cluster labels to a viridis gradient color keyed by sorted label position.
#
# `unique_labels` is any iterable of integer labels. The labels are sorted; the lowest maps to
# one end of the viridis gradient and the highest to the other (position `i/(n-1)`), matching the
# Python pipeline's `_CLUSTER_CMAP = "viridis"`. A single label maps to the gradient midpoint. A
# `-1` noise label, if present, stays light gray. Returns a `Dict` from label to color. Pure and
# deterministic given the labels.
function _cluster_color_map(unique_labels)
    uniq = sort(unique(Int.(collect(unique_labels))))
    grad = [c for c in uniq if c != -1]
    n = length(grad)
    cmap = plt.cgrad(:viridis)
    color_map = Dict{Int,typeof(_NOISE_COLOR)}()
    (-1 in uniq) && (color_map[-1] = _NOISE_COLOR)
    for (i, lbl) in enumerate(grad)
        t = n <= 1 ? 0.5 : (i - 1) / (n - 1)
        color_map[lbl] = cmap[t]
    end
    return color_map
end

# A categorical palette of `k` evenly-spaced, equally-bright hues, approximating seaborn's
# "husl" palette used for the mini-cluster nodes. Returns a vector of `k` RGB colors.
function _husl_palette(k::Integer)
    return [convert(plt.RGB, plt.HSL(360 * (i - 1) / max(k, 1), 0.9, 0.65)) for i = 1:k]
end

# Minimum-spanning-tree edges of a symmetric `K×K` distance matrix (Prim's algorithm).
#
# Returns a vector of `(i, j, weight)` tuples (`i < j`), the `K-1` lightest edges that keep every
# node connected - the structural backbone used by the mini-cluster graph. Pure stdlib; no
# dependency beyond the distance matrix. An empty or singleton matrix yields no edges.
function _mst_edges(dist::AbstractMatrix)
    k = size(dist, 1)
    k < 2 && return Tuple{Int,Int,Float64}[]
    in_tree = falses(k)
    in_tree[1] = true
    edges = Tuple{Int,Int,Float64}[]
    for _ = 1:(k-1)
        best = (0, 0, Inf)
        for i = 1:k
            in_tree[i] || continue
            for j = 1:k
                in_tree[j] && continue
                if dist[i, j] < best[3]
                    best = (i, j, Float64(dist[i, j]))
                end
            end
        end
        best[1] == 0 && break
        i, j, w = best
        in_tree[j] = true
        push!(edges, (min(i, j), max(i, j), w))
    end
    return edges
end

# Average-linkage RMSD between sub-clusters of one KMeans cluster.
#
# `dist` is the member `N×N` pairwise RMSD matrix and `sub_labels` the per-member sub-cluster
# label. Returns `(uniq, mean)` where `uniq` are the sorted unique sub-labels and `mean` the
# `K×K` matrix of average between-sub-cluster RMSD (the diagonal is the mean within-cluster RMSD).
function _intercluster_rmsd(dist::AbstractMatrix, sub_labels::AbstractVector)
    uniq = sort(unique(sub_labels))
    k = length(uniq)
    mean = zeros(Float64, k, k)
    for a = 1:k, b = 1:k
        ia = findall(==(uniq[a]), sub_labels)
        ib = findall(==(uniq[b]), sub_labels)
        mean[a, b] = Statistics.mean(dist[ia, ib])
    end
    return uniq, mean
end

# Cluster-overview scatter: the PCA projection colored by KMeans cluster.
#
# Points are colored by their cluster label as a viridis gradient with a right-hand colorbar (the
# Slice 8 styling), not a categorical legend.
function _plot_cluster_overview(projection::AbstractMatrix, labels::AbstractVector, path)
    p = plt.scatter(
        projection[1, :],
        projection[2, :];
        marker_z = labels,
        color = :viridis,
        colorbar = true,
        colorbar_title = "KMeans cluster",
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
        xlabel = "PC 1",
        ylabel = "PC 2",
        title = "Cluster overview (PCA)",
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.savefig(p, path)
    return path
end

# Adaptive reference scatter: apo-vs-holo (2 refs), embedding-vs-reference (1 ref) or the
# plain MDS embedding (0 refs), colored by KMeans cluster as a viridis gradient with a right-hand
# colorbar. `embedding` is the `2×N` MDS layout used in the 0/1-reference cases; `rmsd_apo`/
# `rmsd_holo` carry the per-conformer reference RMSDs.
function _plot_reference_scatter(
    ref_mode::Integer,
    labels::AbstractVector,
    embedding,
    rmsd_apo,
    rmsd_holo,
    path,
)
    common = (
        marker_z = labels,
        color = :viridis,
        colorbar = true,
        colorbar_title = "KMeans cluster",
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    if ref_mode == 2
        p = plt.scatter(
            rmsd_holo,
            rmsd_apo;
            xlabel = "RMSD to holo reference (Å)",
            ylabel = "RMSD to apo reference (Å)",
            title = "Apo vs holo RMSD",
            common...,
        )
    elseif ref_mode == 1
        p = plt.scatter(
            embedding[1, :],
            rmsd_apo;
            xlabel = "Embedding dimension 1",
            ylabel = "RMSD to apo reference (Å)",
            title = "Embedding vs reference RMSD",
            common...,
        )
    else
        p = plt.scatter(
            embedding[1, :],
            embedding[2, :];
            xlabel = "MDS dimension 1",
            ylabel = "MDS dimension 2",
            title = "Conformer embedding",
            common...,
        )
    end
    plt.savefig(p, path)
    return path
end

# Score-colored RMSD scatter: the PCA projection of the scored conformers colored by score.
function _plot_rmsd_by_score(projection::AbstractMatrix, scores::AbstractVector, path)
    p = plt.scatter(
        projection[1, :],
        projection[2, :];
        marker_z = scores,
        markersize = 5,
        markerstrokewidth = 0,
        colorbar_title = "score",
        xlabel = "PC 1",
        ylabel = "PC 2",
        title = "Embedding colored by score",
        legend = false,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.savefig(p, path)
    return path
end

# KMeans clustering with DeepAccNet filtering: the cluster-overview scatter (viridis gradient)
# showing the filter result - conformers in surviving clusters keep their cluster color, those in
# filtered-out clusters are de-emphasized in light gray. `surviving` is the set of surviving
# cluster labels. Written under `deepaccnet/`. Mirrors the cluster-overview styling (Slice 8 B).
function _plot_kmeans_deepaccnet_filtered(
    projection::AbstractMatrix,
    labels::AbstractVector,
    surviving,
    path,
)
    keep = Set(surviving)
    color_map = _cluster_color_map(unique(labels))
    point_colors = [l in keep ? color_map[l] : _NOISE_COLOR for l in labels]
    n_surv = count(in(keep), labels)
    p = plt.scatter(
        projection[1, :],
        projection[2, :];
        color = point_colors,
        markersize = 5,
        markerstrokewidth = 0,
        markeralpha = 0.7,
        label = "",
        xlabel = "PC 1",
        ylabel = "PC 2",
        title = "KMeans clusters after DeepAccNet filter ($n_surv surviving)",
        legend = false,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.savefig(p, path)
    return path
end

# Dendrogram of one cluster's agglomerative tree with a dashed line at the `threshold` cut.
#
# `tree` is a `Clustering.Hclust`; its `merges` encode children (a negative entry is a leaf, a
# positive entry a previous merge row) and `heights` the merge heights. A `nothing` tree (a
# cluster too small to merge) draws a single leaf marker so the figure still appears.
function _plot_dendrogram(tree, threshold::Real, title::AbstractString, path)
    if tree === nothing
        p = plt.scatter(
            [1.0],
            [0.0];
            legend = false,
            xlabel = "conformer (leaf order)",
            ylabel = "merge height (Å)",
            title = title,
            size = _FIG_SIZE,
            dpi = _FIG_DPI,
            left_margin = _FIG_MARGIN,
            right_margin = _FIG_MARGIN,
            top_margin = _FIG_MARGIN,
            bottom_margin = _FIG_MARGIN,
        )
        plt.savefig(p, path)
        return path
    end

    merges = tree.merges
    heights = tree.heights
    order = tree.order
    n = length(order)
    xpos = zeros(Float64, n)
    for (i, leaf) in enumerate(order)
        xpos[leaf] = i
    end
    merge_x = zeros(Float64, n - 1)
    node_xy(child) = child < 0 ? (xpos[-child], 0.0) : (merge_x[child], heights[child])

    xs = Float64[]
    ys = Float64[]
    for m = 1:(n-1)
        x1, y1 = node_xy(merges[m, 1])
        x2, y2 = node_xy(merges[m, 2])
        h = heights[m]
        merge_x[m] = (x1 + x2) / 2
        append!(xs, [x1, x1, x2, x2, NaN])
        append!(ys, [y1, h, h, y2, NaN])
    end

    p = plt.plot(
        xs,
        ys;
        color = :steelblue,
        legend = false,
        xlabel = "conformer (leaf order)",
        ylabel = "merge height (Å)",
        title = title,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.hline!(p, [threshold]; color = :red, linestyle = :dash)
    plt.savefig(p, path)
    return path
end

# Adaptive sub-cluster scatter: the MDS embedding of one cluster's conformers colored by their
# 1-based sub-cluster label as a viridis gradient with a right-hand colorbar.
function _plot_subcluster_scatter(
    embedding::AbstractMatrix,
    sub_labels::AbstractVector,
    title::AbstractString,
    path,
)
    p = plt.scatter(
        embedding[1, :],
        embedding[2, :];
        marker_z = sub_labels,
        color = :viridis,
        colorbar = true,
        colorbar_title = "sub-cluster",
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
        xlabel = "MDS dimension 1",
        ylabel = "MDS dimension 2",
        title = title,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.savefig(p, path)
    return path
end

# Mini-cluster graph for one KMeans cluster: each node is an agglomerative sub-cluster, laid out
# by 2-D classical MDS of the K×K average-linkage RMSD matrix so screen distance approximates
# inter-cluster RMSD; edges are the minimum-spanning-tree backbone, colored and weighted by the
# between-sub-cluster RMSD (closer = darker/thicker) with a matching colorbar. Skipped (returns
# `nothing`) when the cluster has fewer than two sub-clusters. Mirrors the Python
# `_plot_minicluster_graph` (the node layout uses classical MDS rather than the Python iterative
# MDS, keeping the dimensionality reduction on stdlib linear algebra).
function _plot_minicluster_graph(
    dist::AbstractMatrix,
    sub_labels::AbstractVector,
    title::AbstractString,
    path,
)
    uniq, inter = _intercluster_rmsd(dist, sub_labels)
    k = length(uniq)
    k < 2 && return nothing

    symmetric = (inter .+ inter') ./ 2
    for i = 1:k
        symmetric[i, i] = 0.0
    end
    symmetric[symmetric .< 0] .= 0.0
    coords = _classical_mds(symmetric, 2)
    edges = _mst_edges(symmetric)
    sizes = [count(==(u), sub_labels) for u in uniq]
    palette = _husl_palette(k)        # categorical node colors, matching seaborn "husl"
    edge_cmap = plt.cgrad(:viridis; rev = true)   # viridis_r, matching the Python edge colormap

    weights = isempty(edges) ? [1.0] : [e[3] for e in edges]
    wmin, wmax = minimum(weights), max(maximum(weights), minimum(weights) + 1e-6)

    p = plt.plot(;
        legend = :outerbottom,
        legendtitle = "mini-cluster (N = #structures)",
        legend_column = min(k, 8, max(3, (k + 4) ÷ 5)),
        xlabel = "MDS-1 (Å)",
        ylabel = "MDS-2 (Å)",
        title = "$title\n$k mini-clusters; edge color = avg-linkage RMSD (see bar); " *
                "MST backbone ($(length(edges)) edges)",
        aspect_ratio = :equal,
        grid = true,
        gridstyle = :dash,
        gridalpha = 0.3,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )

    # Edges: viridis_r colored and weighted by between-sub-cluster RMSD (closer = darker/thicker).
    for (i, j, w) in edges
        t = (w - wmin) / (wmax - wmin)
        plt.plot!(
            p,
            [coords[1, i], coords[1, j]],
            [coords[2, i], coords[2, j]];
            color = edge_cmap[t],
            linewidth = 0.6 + 2.2 * (1 - t),
            alpha = 0.85,
            label = "",
        )
    end
    # A hidden 2-point series carries the edge-RMSD range so a colorbar is drawn.
    if !isempty(edges)
        plt.scatter!(
            p,
            [coords[1, 1], coords[1, 1]],
            [coords[2, 1], coords[2, 1]];
            marker_z = [wmin, wmax],
            color = edge_cmap,
            markeralpha = 0,
            colorbar = true,
            colorbar_title = "edge = avg-linkage RMSD between mini-clusters (Å)",
            label = "",
        )
    end
    # Nodes: one series per mini-cluster so the legend lists each with its member count.
    for idx = 1:k
        plt.scatter!(
            p,
            [coords[1, idx]],
            [coords[2, idx]];
            color = palette[idx],
            markersize = 10,
            markerstrokewidth = 0.8,
            markerstrokecolor = :black,
            label = "$(uniq[idx]) (N=$(sizes[idx]))",
            series_annotations = plt.text(string(uniq[idx]), 8, :black),
        )
    end
    plt.savefig(p, path)
    return path
end

# Query-path view (cluster-colored): rank every KMeans cluster along x by centroid-to-centroid
# CA-RMSD to the query cluster (the query cluster is rank 0 at y = 0), with y the RMSD to the
# query cluster. Nodes are colored by the KMeans cluster gradient; reference-assigned clusters
# are overlaid with a star. Mirrors the Python `build_query_path_cluster_colored_plot`.
function _plot_query_path(
    inter_rmsd::AbstractMatrix,
    query_cluster::Integer,
    cluster_ids::AbstractVector,
    ref_clusters,
    path,
)
    n = length(cluster_ids)
    qpos = findfirst(==(query_cluster), cluster_ids)
    others = sort(
        [c for c in cluster_ids if c != query_cluster];
        by = c -> inter_rmsd[qpos, findfirst(==(c), cluster_ids)],
    )
    visited = vcat([query_cluster], others)
    xs = collect(0:(n-1))
    ys = [inter_rmsd[qpos, findfirst(==(c), cluster_ids)] for c in visited]

    colors = _cluster_color_map(cluster_ids)
    p = plt.plot(
        xs,
        ys;
        color = :gray,
        linewidth = 1.2,
        alpha = 0.7,
        legend = false,
        xlabel = "cluster rank by RMSD to query (0 = query cluster)",
        ylabel = "CA-RMSD to query cluster (Å)",
        title = "Query distance profile (cluster colors)",
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    for (i, c) in enumerate(visited)
        plt.scatter!(
            p,
            [xs[i]],
            [ys[i]];
            color = colors[c],
            markersize = c == query_cluster ? 11 : 8,
            markerstrokewidth = c == query_cluster ? 2.5 : 0.8,
            markerstrokecolor = :black,
            series_annotations = plt.text(string(c), 7, :black),
        )
    end
    if ref_clusters !== nothing
        for (name, cid) in ref_clusters
            (cid === nothing || !(cid in cluster_ids)) && continue
            i = findfirst(==(cid), visited)
            plt.scatter!(
                p,
                [xs[i]],
                [ys[i]];
                color = colors[cid],
                markershape = :star5,
                markersize = 16,
                markerstrokewidth = 1.8,
                markerstrokecolor = :black,
            )
        end
    end
    plt.savefig(p, path)
    return path
end

# Render the kept figure set for one run into the per-cluster output tree, adapting to the
# reference count and the optional score filter. Returns the vector of written figure paths.
#
# Pure orchestration over the `_plot_*` helpers; computes the PCA/MDS layouts on the fly. Sets
# `GKSwstype` so the GR backend renders offscreen on headless machines. System-level figures land
# at `system_dir`, score-filter figures under `deepaccnet_dir`, and each KMeans cluster's figures
# in its own `kmeans<C>/` folder (carried in `cluster_results`).
function _render_figures(;
    system_dir,
    deepaccnet_dir,
    ref_mode::Integer,
    labels::AbstractVector,
    features::AbstractMatrix,
    rmsd_apo,
    rmsd_holo,
    reference_clusters,
    scored_idx::AbstractVector,
    scores::AbstractVector,
    surviving,
    cluster_results,
    linkage::Symbol,
    threshold::Real,
    query_cluster,
    inter_rmsd,
    cluster_ids,
)
    get!(ENV, "GKSwstype", "100")
    figures = String[]

    # System-level clustering figures.
    projection = _pca(features, 2)
    push!(
        figures,
        _plot_cluster_overview(
            projection,
            labels,
            joinpath(system_dir, "cluster_overview.png"),
        ),
    )

    embedding =
        ref_mode == 2 ? zeros(2, length(labels)) :
        _classical_mds(_feature_distance_matrix(features), 2)
    push!(
        figures,
        _plot_reference_scatter(
            ref_mode,
            labels,
            embedding,
            rmsd_apo,
            rmsd_holo,
            joinpath(system_dir, "reference_scatter.png"),
        ),
    )

    # Query-path view - only when a query cluster was resolved (C4).
    if query_cluster !== nothing
        ref_map =
            reference_clusters === nothing ? nothing :
            Dict(
                string(r) => c for
                (r, c) in zip(reference_clusters.Reference, reference_clusters.Cluster)
            )
        push!(
            figures,
            _plot_query_path(
                inter_rmsd,
                query_cluster,
                cluster_ids,
                ref_map,
                joinpath(system_dir, "query_path.png"),
            ),
        )
    end

    # Score-filter figures (only when the filter ran).
    if deepaccnet_dir !== nothing
        push!(
            figures,
            _plot_rmsd_by_score(
                projection[:, scored_idx],
                scores,
                joinpath(deepaccnet_dir, "rmsd_scatter_by_score.png"),
            ),
        )
        push!(
            figures,
            _plot_kmeans_deepaccnet_filtered(
                projection,
                labels,
                surviving,
                joinpath(deepaccnet_dir, "kmeans_deepaccnet_filtered.png"),
            ),
        )
    end

    # Per-cluster agglomerative figures, one set per `kmeans<C>/` folder.
    for cr in cluster_results
        dist = cr.dist
        sub_labels = cr.sub_labels
        cluster_dir = cr.cluster_dir
        title = "Cluster $(cr.lbl) @ $(threshold) Å cut"

        cluster_embedding = _classical_mds(dist, 2)
        push!(
            figures,
            _plot_subcluster_scatter(
                cluster_embedding,
                sub_labels,
                title,
                joinpath(cluster_dir, "subcluster_scatter.png"),
            ),
        )

        tree = size(dist, 1) >= 2 ? Clustering.hclust(dist; linkage = linkage) : nothing
        push!(
            figures,
            _plot_dendrogram(
                tree,
                threshold,
                title,
                joinpath(cluster_dir, "dendrogram.png"),
            ),
        )
        graph = _plot_minicluster_graph(
            dist,
            sub_labels,
            "Cluster $(cr.lbl) mini-clusters",
            joinpath(cluster_dir, "minicluster_graph.png"),
        )
        graph === nothing || push!(figures, graph)
    end

    return figures
end

# Centroid-to-centroid CA-RMSD between every occupied KMeans cluster.
#
# `features` is the `(3L)×N` aligned-coordinate matrix and `labels` the 0-based cluster labels.
# Returns `(cluster_ids, inter_rmsd)` where `cluster_ids` are the sorted occupied labels and
# `inter_rmsd` the matching `C×C` matrix of centroid euclidean distance divided by `sqrt(n_atoms)`
# (the CA-RMSD between cluster centroids), mirroring the Python query-path metric.
function _cluster_centroid_rmsd(
    features::AbstractMatrix,
    labels::AbstractVector,
    n_atoms::Integer,
)
    cluster_ids = sort(unique(labels))
    c = length(cluster_ids)
    centroids = Matrix{Float64}(undef, size(features, 1), c)
    for (j, lbl) in enumerate(cluster_ids)
        idx = findall(==(lbl), labels)
        centroids[:, j] = vec(Statistics.mean(features[:, idx]; dims = 2))
    end
    inter = zeros(Float64, c, c)
    for i = 1:c, j = (i+1):c
        d = sqrt(sum(abs2, view(centroids, :, i) .- view(centroids, :, j))) / sqrt(n_atoms)
        inter[i, j] = d
        inter[j, i] = d
    end
    return cluster_ids, inter
end

# Resolve the query input to a KMeans cluster id.
#
# `query` is either a conformer name present in `names` (its cluster is returned directly) or a
# path to a structure file. For a path, the query is assigned to the cluster whose conformers are
# closest to it on average, using the same sequence-aware RMSD as the reference assignment (which
# tolerates a query whose residue numbering differs from the ensemble). Returns the cluster id.
# Raises a clear error when a name is not found and the value is not an existing file.
function _resolve_query_cluster(
    query::AbstractString,
    names::AbstractVector,
    labels::AbstractVector,
    conformers::AbstractVector,
)
    pos = findfirst(==(query), names)
    pos !== nothing && return labels[pos]
    isfile(query) ||
        error("query '$query' is neither a conformer name nor an existing structure file")

    query_residues = _read_ca_residues(query)
    conformer_residues = [_read_ca_residues(path) for (_name, path) in conformers]
    rmsds = _reference_rmsds(query_residues, conformer_residues)
    cluster, _mean = _assign_reference_cluster(rmsds, labels)
    return cluster
end

# Public entry point
# -----

"""
    cluster_conformers(system, apo_ref="", holo_ref="";
        refs_dir="refs", data_root="/data/alphaconformers/pdb", out_dir="results",
        kmeans_k=30, threshold=1.0, linkage=:average,
        score_table=nothing, score_col="cb-lddt", surviving_rel=25, surviving_frac=75,
        seed=42, make_plots=true, query=nothing)

Cluster an AlphaConformer prediction ensemble by shape, write the base KMeans tables, and
(when references are supplied) report each reference against the cluster it falls into.

It discovers the system's conformers, keeps the CA positions shared by all of them,
superimposes the coordinates, runs a seeded KMeans at `kmeans_k`, and writes the clustering
table and the intra-cluster RMSD table at the per-system output root. The reference mode is
auto-detected from the stems supplied: 0 references (reference-free), 1 (apo only), or 2 (apo
and holo). With one or two references it also computes the sequence-aware RMSD from every
conformer to each reference, assigns each reference to the cluster whose conformers are closest
to it on average, and writes the apo/holo RMSD table and the reference cluster-assignment table.
When a DeepAccNet score table is supplied it runs the optional score filter: for each
relative/fraction cut it finds the surviving clusters (a cluster survives when at most the
fraction of its members fall below the relative score cutoff) and writes a surviving-cluster
table per cut under a `deepaccnet/` folder; with no score table this stage is skipped. It then
always sub-clusters: within each KMeans cluster it builds the pairwise RMSD matrix of the
members, cuts an agglomerative tree at `threshold`, and writes that cluster's `members.csv` into
its own `kmeans<C>/` folder, while a flat top-level `agglomerative_assignments.csv` keeps the
single-table view. When the score filter ran, only the surviving KMeans clusters get a
`kmeans<C>/` folder. Unless `make_plots` is disabled, it also writes the figures for each level,
adapting to the reference count (0/1/2). The query-path view is added whenever a reference is
given (it tracks the apo reference) or an explicit `query` is supplied.

# Arguments

- `system::String`: System identifier; names the input subdirectory under `data_root` and the
  output subdirectory under `out_dir`. Must not be blank after whitespace stripping.
- `apo_ref::String`: (Optional) apo reference stem, e.g. `"1ABC_A"`. An empty/blank stem
  selects reference-free or, with `holo_ref`, is rejected. Defaults to `""`.
- `holo_ref::String`: (Optional) holo reference stem. If given, `apo_ref` must also be given.
  Defaults to `""`.

# Keywords

- `refs_dir::String`: Directory holding the apo/holo reference files. Each stem is resolved
  against this directory. Defaults to `"refs"`.
- `data_root::String`: Root directory holding one AlphaConformers output tree per `<system>`
  subfolder; predictions are read from the `models/` folders inside it. Defaults to
  `"/data/alphaconformers/pdb"`, matching where the package writes prediction outputs.
- `out_dir::String`: Root directory for output. Defaults to `"results"`.
- `kmeans_k::Int`: Number of clusters. Must be positive and at most the conformer count.
  Defaults to `30`.
- `threshold::Float64`: Distance threshold (Å) at which the agglomerative tree is cut into
  sub-clusters. Must be positive. Defaults to `1.0`.
- `linkage::Symbol`: Linkage method for the agglomerative sub-clustering (e.g. `:average`,
  `:single`, `:complete`). Defaults to `:average`.
- `score_table`: Optional DeepAccNet score table, either a `DataFrame` or a path to a CSV. When
  supplied, the score filter runs; when `nothing` (the default) the filter is skipped. Must hold
  a `Name` column matching the conformers and the `score_col` score column. Conformers absent
  from the table are dropped from the filter.
- `score_col::String`: Name of the score column read from `score_table`. Defaults to `"cb-lddt"`.
- `surviving_rel::Int`: Relative percentile cut (`0..100`) defining the score cutoff for the
  surviving-cluster set that feeds the agglomerative stage. Only this configured cut is computed.
  Defaults to `25`.
- `surviving_frac::Int`: Fraction-removed cut in percent (`0..100`); a cluster survives when at
  most this fraction of its members fall below the relative cutoff. Defaults to `75`.
- `seed::Int`: Random seed for the KMeans RNG, ensuring reproducible clustering. Defaults to
  `42`.
- `make_plots::Bool`: Whether to write the figures alongside the tables. Defaults to `true`.
- `query`: Optional query input for the query-path view - either a conformer name (as it appears
  in the clustering table) or a path to a structure file. When `nothing` (the default) the query
  falls back to the apo reference if any reference was given: with one reference that reference is
  used, with two references the apo is used. So a run with references always writes the query-path
  figure tracking apo, and passing `query` explicitly overrides this. In reference-free mode with
  no `query` the query-path figure is skipped.

# Returns

A named tuple `(; clustering, cluster_rmsd, system_dir, reference_rmsd, reference_clusters,
deepaccnet_dir, surviving, hierarchical, cluster_dirs, figures)`:

- `clustering`: `DataFrame` with columns `Name` and `KMeans_K<k>` (labels in `0:k-1`).
- `cluster_rmsd`: `DataFrame` with columns `Method`, `Cluster`, `Size`, `Mean_RMSD_Angstrom`.
- `system_dir`: path to the per-system output root holding the flat tables, the system-level
  figures, the `kmeans<C>/` folders and the optional `deepaccnet/` folder.
- `reference_rmsd`: `nothing` in reference-free mode; otherwise a `DataFrame` with a `Name`
  column and a `RMSD_Apo` and/or `RMSD_Holo` column, one row per conformer.
- `reference_clusters`: `nothing` in reference-free mode; otherwise a `DataFrame` with columns
  `Reference`, `Stem`, `Cluster`, `Mean_RMSD_Angstrom`, one row per supplied reference.
- `deepaccnet_dir`: `nothing` when no score table is supplied; otherwise the path to the written
  `deepaccnet/` directory.
- `surviving`: `nothing` when no score table is supplied; otherwise the sorted vector of
  surviving cluster labels for the configured `surviving_rel`/`surviving_frac` cut.
- `hierarchical`: flat `DataFrame` with columns `Name`, `KMeans_K<k>`, `Sub_Cluster` (1-based
  label within the KMeans cluster) and `Mini_Cluster` (the `c<label>_s<sub>` mini-cluster name).
  One row per conformer, or, when the score filter ran, one row per conformer in a surviving
  cluster.
- `cluster_dirs`: vector of the written `kmeans<C>/` folder paths, one per occupied (or
  surviving) KMeans cluster.
- `figures`: vector of paths to the figures written; empty when `make_plots` is `false`.

# Throws

- `ErrorException`: if `system` is blank, if `kmeans_k <= 0` or `threshold <= 0`, if
  `holo_ref` is given without `apo_ref`, if a reference stem cannot be resolved under
  `refs_dir`, if no conformers are found, if the conformers share no common CA
  position, if `kmeans_k` exceeds the conformer count, or, when a score table is supplied, if it
  is missing the `Name`/`score_col` column, if it matches no conformer, if `surviving_rel` or
  `surviving_frac` is outside `0..100`, or if a `query` is given that is neither a known
  conformer name nor an existing structure file.

# Behavior

1. Validate inputs and detect the reference mode (0 = none, 1 = apo, 2 = apo and holo).
2. Discover the prediction structures under `data_root/<system>` - only the files in `models/`
   folders (the AlphaConformers prediction location), so input templates are skipped - parse
   their CA traces, and keep the CA positions shared by all conformers.
3. Superimpose every conformer onto the first and run seeded KMeans at `kmeans_k`.
4. Compute the mean intra-cluster RMSD to each cluster centroid and write
   `aligned_clustering_results.csv` and `aligned_cluster_rmsd.csv` at `out_dir/<system>/`.
5. With one or two references, resolve each stem against `refs_dir`, compute the sequence-aware
   RMSD from every conformer to each reference, assign each reference to the cluster with the
   lowest mean RMSD, and write `dist_external_rmsds.csv` and `reference_clusters.csv` at
   `out_dir/<system>/`.
6. When a score table is supplied, align its scores to the conformers (dropping unscored ones),
   find the surviving clusters for the configured `surviving_rel`/`surviving_frac` cut, and write
   `surviving_rel<rel>_frac<frac>.csv` under `out_dir/<system>/deepaccnet/surviving/`. Skipped
   when no score table is given.
7. For each occupied (or, when the filter ran, surviving) KMeans cluster `C`, build the pairwise
   RMSD matrix of its members, cut the agglomerative tree at `threshold` with the chosen
   `linkage`, and write that cluster's `members.csv` into `out_dir/<system>/kmeans<C>/`. A flat
   top-level `agglomerative_assignments.csv` collects all sub-cluster assignments.
8. When `make_plots` is set, write the figures at each level: the
   cluster-overview scatter, the adaptive reference scatter, and the query-path view (whenever a
   reference is given, or an explicit `query`) at the system root; the score-colored embedding and
   DeepAccNet-filtered
   scatter under `deepaccnet/` when the filter ran; and the dendrogram, sub-cluster embedding and
   mini-cluster graph in each `kmeans<C>/`.
9. Return the tables, the output directories and the written figure paths.
"""
function cluster_conformers(
    system::String,
    apo_ref::String = "",
    holo_ref::String = "";
    refs_dir::String = "refs",
    data_root::String = "/data/alphaconformers/pdb",
    out_dir::String = "results",
    kmeans_k::Int = 30,
    threshold::Float64 = 1.0,
    linkage::Symbol = :average,
    score_table = nothing,
    score_col::String = "cb-lddt",
    surviving_rel::Int = 25,
    surviving_frac::Int = 75,
    seed::Int = 42,
    make_plots::Bool = true,
    query::Union{Nothing,AbstractString} = nothing,
)
    # Validate and normalize inputs.
    stripped_system = strip(system)
    isempty(stripped_system) && error("system must not be blank")

    kmeans_k > 0 || error("kmeans_k must be positive")
    threshold > 0 || error("threshold must be positive")

    apo_normalized = _normalize_reference(apo_ref)
    holo_normalized = _normalize_reference(holo_ref)

    # Detect reference mode (0 = none, 1 = apo, 2 = apo and holo).
    ref_mode = _reference_mode(apo_normalized, holo_normalized)

    # Discover and parse the ensemble's conformers (CA coordinates).
    conformers = _discover_conformers(data_root, String(stripped_system))
    names = first.(conformers)
    parsed = [_parse_ca(path) for (_name, path) in conformers]
    resids_list = first.(parsed)

    if kmeans_k > length(conformers)
        error(
            "kmeans_k ($kmeans_k) exceeds the number of conformers " *
            "($(length(conformers))); choose a smaller kmeans_k.",
        )
    end

    # Restrict to the common CA positions and superimpose.
    common = _common_ca_positions(resids_list)
    restricted =
        [_restrict_to_common(resids, coords, common) for (resids, coords) in parsed]
    aligned = _superimpose_to_first(restricted)
    features = _feature_matrix(aligned)

    # Seeded KMeans and intra-cluster RMSD.
    labels = _kmeans_cluster(features, kmeans_k; seed = seed)
    rmsd = _intra_cluster_rmsd(features, labels, length(common))

    # Per-system output root: the top-level flat tables and the system-level figures live here,
    # alongside the per-cluster `kmeans<C>/` folders and the optional `deepaccnet/` folder.
    system_dir = joinpath(out_dir, String(stripped_system))
    clustering, cluster_rmsd =
        _write_kmeans_tables(system_dir, names, labels, kmeans_k, rmsd)

    # Reference handling (apo/holo): sequence-aware RMSD to every conformer and the cluster
    # each reference falls into (lowest mean RMSD). Skipped entirely in reference-free mode.
    reference_rmsd = nothing
    reference_clusters = nothing
    rmsd_apo = Float64[]
    rmsd_holo = Float64[]
    if ref_mode > 0
        conformer_residues = [_read_ca_residues(path) for (_name, path) in conformers]

        refs = Tuple{String,String,Vector{Float64}}[]
        rmsd_apo = _reference_rmsds(
            _read_ca_residues(_resolve_reference(refs_dir, apo_normalized)),
            conformer_residues,
        )
        push!(refs, ("apo", apo_normalized, rmsd_apo))
        if ref_mode == 2
            rmsd_holo = _reference_rmsds(
                _read_ca_residues(_resolve_reference(refs_dir, holo_normalized)),
                conformer_residues,
            )
            push!(refs, ("holo", holo_normalized, rmsd_holo))
        end

        reference_rmsd, reference_clusters =
            _write_reference_tables(system_dir, names, labels, refs)
    end

    # DeepAccNet score filter (Stage 2). Runs only when a score table is supplied. It computes the
    # surviving clusters for the configured `surviving_rel`/`surviving_frac` cut (a cluster
    # survives when at most the fraction of its members fall below the relative cutoff) and writes
    # that single surviving-cluster table under `deepaccnet/`; the cut feeds the agglomerative
    # stage below. Skipped cleanly when no score table is given. Reference-free: it needs scores
    # only.
    deepaccnet_dir = nothing
    surviving = nothing
    scored_idx = Int[]
    scores = Float64[]
    if score_table !== nothing
        0 <= surviving_rel <= 100 ||
            error("surviving_rel must be a percentile in 0..100, got $surviving_rel")
        0 <= surviving_frac <= 100 ||
            error("surviving_frac must be a percentage in 0..100, got $surviving_frac")
        deepaccnet_dir = joinpath(system_dir, "deepaccnet")
        scored_idx, scores = _align_scores(score_table, names; score_col = score_col)
        surviving = _write_surviving_table(
            deepaccnet_dir,
            names[scored_idx],
            labels[scored_idx],
            scores,
            kmeans_k,
            surviving_rel,
            surviving_frac,
        )
    end

    # Agglomerative sub-clustering (always runs). Within each KMeans cluster, build the pairwise
    # RMSD matrix of its members, cut the agglomerative tree at `threshold`, and record a
    # per-conformer sub-cluster label. When the score filter ran, only the scored conformers in
    # the configured surviving clusters are sub-clustered; otherwise the full KMeans clustering.
    # Each cluster's outputs land in its own `kmeans<C>/` folder, with a leaf folder per
    # sub-cluster; a flat top-level `agglomerative_assignments.csv` keeps the single-table view.
    if score_table === nothing
        work_idx = collect(eachindex(names))
    else
        keep = Set(surviving)
        work_idx = [i for i in scored_idx if labels[i] in keep]
    end
    work_names = names[work_idx]
    work_labels = labels[work_idx]
    work_aligned = aligned[work_idx]

    sub_labels = zeros(Int, length(work_idx))
    cluster_results = NamedTuple[]
    cluster_dirs = String[]
    for lbl in sort(unique(work_labels))
        idx = findall(==(lbl), work_labels)
        dist = _pairwise_rmsd_matrix(work_aligned[idx])
        subs = _agglomerative_subcluster(dist, threshold, linkage)
        sub_labels[idx] .= subs

        cluster_dir = joinpath(system_dir, "kmeans$(lbl)")
        _write_cluster_outputs(
            cluster_dir,
            work_names[idx],
            work_labels[idx],
            subs,
            kmeans_k,
        )
        push!(cluster_dirs, cluster_dir)
        push!(
            cluster_results,
            (; lbl = lbl, dist = dist, sub_labels = subs, cluster_dir = cluster_dir),
        )
    end

    # Flat top-level assignments table (single source for downstream consumers).
    hierarchical = _hierarchical_table(work_names, work_labels, sub_labels, kmeans_k)
    CSV.write(joinpath(system_dir, "agglomerative_assignments.csv"), hierarchical)

    # Query-path inputs (C4). The query defaults to the apo reference when references are given
    # (1 reference: that reference; 2 references: the apo), so the query-path view tracks apo
    # without an extra argument. An explicit `query` keyword overrides this default. Centroid
    # distances are computed over the full KMeans clustering so the ranking covers every occupied
    # cluster. For the defaulted apo query the cluster is apo's already-computed assigned cluster
    # (row 1 of `reference_clusters`), avoiding a second RMSD pass.
    query_cluster = nothing
    cluster_ids = Int[]
    inter_rmsd = zeros(0, 0)
    if query !== nothing
        cluster_ids, inter_rmsd = _cluster_centroid_rmsd(features, labels, length(common))
        query_cluster = _resolve_query_cluster(query, names, labels, conformers)
    elseif ref_mode > 0
        cluster_ids, inter_rmsd = _cluster_centroid_rmsd(features, labels, length(common))
        query_cluster = reference_clusters.Cluster[1]
    end

    # Plotting (Stage: figures). Kept out of the algorithmic functions above; adapts to the
    # reference count and the optional score filter. Disabled with `make_plots = false`.
    figures = String[]
    if make_plots
        figures = _render_figures(;
            system_dir,
            deepaccnet_dir,
            ref_mode,
            labels,
            features,
            rmsd_apo,
            rmsd_holo,
            reference_clusters,
            scored_idx,
            scores,
            surviving,
            cluster_results,
            linkage,
            threshold,
            query_cluster,
            inter_rmsd,
            cluster_ids,
        )
    end

    return (;
        clustering,
        cluster_rmsd,
        system_dir,
        reference_rmsd,
        reference_clusters,
        deepaccnet_dir,
        surviving,
        hierarchical,
        cluster_dirs,
        figures,
    )
end
