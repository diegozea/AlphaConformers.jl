# Conformer clustering pipeline -------------------------------------------------------

# Reference handling
# -----

# Normalize a reference stem: blank/whitespace → nothing, else unchanged.
function _normalize_reference(ref::String)::Union{Nothing,String}
    stripped = strip(ref)
    return isempty(stripped) ? nothing : stripped
end

# Number of references supplied (0, 1 or 2). Either positional may be given on its own; there is
# no required ordering between the two reference slots.
function _reference_mode(ref1::Union{Nothing,String}, ref2::Union{Nothing,String})::Int
    return (ref1 !== nothing) + (ref2 !== nothing)
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
# By default it recurses the system directory and keeps only structure files whose immediate
# parent directory is a `models/` folder (the AlphaConformers prediction location); files under
# template or search-result folders are skipped, so templates are never mistaken for predictions.
# When `subdir` is given it instead globs structure files directly under that relative subfolder
# (which may contain glob wildcards, e.g. `cluster_*/af/predictions/sequences/models`), avoiding a
# full directory walk when the caller already knows the method subfolder. Returns a vector of
# `(name, path)` tuples sorted by `name`, where `name` is the file path relative to the system
# directory (stable across machines) and `path` is the absolute path. Raises a clear error if the
# system directory is missing or holds no prediction structures.
function _discover_conformers(
    data_root::AbstractString,
    system::AbstractString;
    subdir::Union{Nothing,AbstractString} = nothing,
)
    system_dir = joinpath(data_root, system)
    isdir(system_dir) ||
        error("no conformer directory found for system '$system' at '$system_dir'")

    found = Tuple{String,String}[]
    if subdir === nothing
        for (root, _dirs, files) in walkdir(system_dir)
            basename(root) == _MODELS_DIR || continue
            for f in files
                ext = lowercase(splitext(f)[2])
                if ext == ".pdb" || ext == ".cif"
                    path = joinpath(root, f)
                    push!(found, (relpath(path, system_dir), abspath(path)))
                end
            end
        end
    else
        for pattern in ("*.pdb", "*.cif")
            for path in glob(joinpath(subdir, pattern), system_dir)
                isfile(path) || continue
                push!(found, (relpath(path, system_dir), abspath(path)))
            end
        end
    end

    isempty(found) && error(
        subdir === nothing ?
        "no .pdb/.cif prediction structures found for system '$system' under '$system_dir'; " *
        "predictions are read only from `$_MODELS_DIR/` folders (the AlphaConformers " *
        "prediction location), so input templates are not treated as conformers." :
        "no .pdb/.cif structures found for system '$system' under '$system_dir/$subdir'",
    )
    sort!(found; by = first)
    return found
end

# Common CA selection (sequence alignment, ensemble-wide)
# -----

# Build one consistent common-Cα set across the whole ensemble by sequence alignment.
#
# `residues_list` is a vector of CA-only `PDBResidue` vectors, one per conformer. Every conformer
# is aligned to the first via `structural_alignment` (the single source of truth for structural
# distance, which matches residues by sequence rather than by residue number), and only the
# reference residue positions matched in *every* conformer are kept, in ascending reference order.
# Returns a vector with each conformer's residues restricted to that consistent set (equal length
# across all conformers). This guards against unreliable residue numbers and against the per-pair
# common sets differing across pairs. Raises when a conformer cannot be aligned or the conformers
# share no common residue.
function _ensemble_common_residues(residues_list::AbstractVector)
    isempty(residues_list) && error("no conformers to find common CA residues for")
    reference = residues_list[1]
    # Per-conformer map from a reference residue index to this conformer's residue index.
    maps = Vector{Dict{Int,Int}}(undef, length(residues_list))
    common = Set{Int}(eachindex(reference))
    for (c, residues) in enumerate(residues_list)
        if c == 1
            maps[1] = Dict(i => i for i in eachindex(reference))
        else
            result = structural_alignment(reference, residues)
            result === nothing &&
                error("could not align conformer $c to the reference; cannot cluster")
            maps[c] = Dict(i_ref => i_conf for (i_ref, i_conf) in result[3])
            intersect!(common, keys(maps[c]))
        end
    end
    common_sorted = sort(collect(common))
    isempty(common_sorted) &&
        error("conformers share no common CA residue positions; cannot cluster")
    return [
        [residues_list[c][maps[c][i]] for i in common_sorted] for
        c in eachindex(residues_list)
    ]
end

# Superimpose every conformer's common-residue set onto the first and return the coordinates.
#
# `common_list` is a vector of equal-length CA-only `PDBResidue` vectors (the ensemble-common set).
# Returns a matching vector of `L×3` coordinate matrices in the first conformer's centered frame,
# using the package's MIToS superposition and `CAmatrix`.
function _superimpose_common(common_list::AbstractVector)
    reference = common_list[1]
    ref_centered, _, _ = MIToS.PDB.superimpose(reference, reference)
    aligned = Vector{Matrix{Float64}}(undef, length(common_list))
    aligned[1] = MIToS.PDB.CAmatrix(ref_centered)
    for i = 2:length(common_list)
        _, bsuper, _ = MIToS.PDB.superimpose(reference, common_list[i])
        aligned[i] = MIToS.PDB.CAmatrix(bsuper)
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

# Structure file extensions tried when resolving a reference stem, including the gzip-compressed
# variants (MIToS `read_file` decompresses `.gz` transparently).
const _STRUCTURE_EXTENSIONS = (".pdb", ".cif", ".pdb.gz", ".cif.gz")

# Pick the MIToS structure file format from a path, looking past a trailing `.gz`.
function _structure_format(path::AbstractString)
    name = endswith(lowercase(path), ".gz") ? path[1:(end-3)] : path
    ext = lowercase(splitext(name)[2])
    return ext == ".cif" ? MIToS.PDB.MMCIFFile : MIToS.PDB.PDBFile
end

# Resolve a reference stem (`<PDBID>_<CHAIN>`) to a file under `refs_dir`.
#
# Tries, across the structure extensions (`.pdb`, `.cif`, and their `.gz` variants): the exact stem
# first, then a case-insensitive lowercase-id fallback that lowercases the PDB id but keeps the
# chain letter (`1AKZ_A` → `1akz_A`). Returns the first existing path. Raises a clear error when the
# stem resolves to no file. No chain argument: each reference file is read whole.
function _resolve_reference(refs_dir::AbstractString, stem::AbstractString)
    candidates = String[String(stem)]
    if occursin('_', stem)
        pdb_id, chain_letter = rsplit(String(stem), '_'; limit = 2)
        push!(candidates, "$(lowercase(pdb_id))_$(chain_letter)")
    end
    for cand in candidates, ext in _STRUCTURE_EXTENSIONS
        path = joinpath(refs_dir, cand * ext)
        isfile(path) && return path
    end
    tried = join([cand * ext for cand in candidates for ext in _STRUCTURE_EXTENSIONS], ", ")
    error("could not resolve reference stem '$stem' under '$refs_dir'; tried: $tried")
end

# Read one structure file into its CA-only `PDBResidue`s.
#
# Reads only the alpha carbons up front (`atomname="CA"`; `.cif` → MMCIF, otherwise PDB, looking
# past a trailing `.gz`), then keeps residues that carry a CA atom. Modified residues are left as
# read - `structural_alignment`/`MIToS.PDB.modelled_sequences` map them to their parent when
# building the sequence. No chain argument: the single-chain-per-file convention means the whole
# file is one chain.
function _read_ca_residues(path::AbstractString)
    residues = MIToS.PDB.read_file(path, _structure_format(path); atomname = "CA")
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

# Write the reference RMSD table and the reference cluster-assignment table.
#
# `refs` is a vector of `(label, stem, rmsds)` tuples (`label` is the user-chosen reference name,
# e.g. `"ref1"`/`"ref2"` by default). Writes `dist_external_rmsds.csv` (one row per conformer, a
# `RMSD_<label>` column per reference) and `reference_clusters.csv` (one row per reference), both
# at the per-system output root. Returns both tables.
function _write_reference_tables(
    system_dir::AbstractString,
    names::AbstractVector,
    labels::AbstractVector,
    refs::AbstractVector,
)
    rmsd_table = DataFrames.DataFrame(:Name => collect(names))
    assignment_rows = NamedTuple[]
    for (label, stem, rmsds) in refs
        rmsd_table[!, Symbol("RMSD_$label")] = rmsds
        cluster, mean_rmsd = _assign_reference_cluster(rmsds, labels)
        push!(
            assignment_rows,
            (
                Reference = label,
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
# Returns 1-based labels (`1:k`). The RNG is seeded from `seed` so the same inputs reproduce the
# same clustering.
function _kmeans_cluster(features::AbstractMatrix, k::Integer; seed::Integer)
    rng = Random.MersenneTwister(seed)
    result = Clustering.kmeans(features, k; rng = rng)
    return Int.(Clustering.assignments(result))
end

# Mean intra-cluster RMSD to the cluster centroid, per cluster label.
#
# `features` is the `(3L)×N` aligned-coordinate matrix, `labels` the cluster assignments and
# `n_atoms` the number of CA positions. Singleton clusters get `0.0`.
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
# `aligned_clustering_results.csv` has columns `Name` and `KMeans_K<k>` (labels `1:k`).
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
# `Name` column and the `score_col` score column. The table is matched to `names` with an inner
# join, so conformers absent from the table are dropped and the filter operates only on scored
# conformers. Returns `(scored_idx, scores)`, where `scored_idx` are the positions in `names` that
# have a score (in `names` order) and `scores` the matching score values. Raises an error if a
# required column is missing or if no conformer matches at all.
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

    conformers =
        DataFrames.DataFrame(_position = collect(eachindex(names)), Name = string.(names))
    scores_by_name = DataFrames.DataFrame(
        Name = string.(table.Name),
        _score = Float64.(table[!, score_col]),
    )
    matched = DataFrames.innerjoin(conformers, scores_by_name; on = :Name)
    sort!(matched, :_position)   # keep the `names` order the join does not guarantee

    isempty(matched) &&
        error("score table matched none of the conformers by `Name`; check the score table")
    return matched._position, matched._score
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
# `rel`/`frac` cut and writes the scored conformers in those clusters to `surviving.csv` with
# columns `Name` and `KMeans_K<k>` (matching the kmeans clustering table). Only this single
# configured cut is produced, so the cut values are not encoded in the filename. Returns the sorted
# vector of surviving cluster labels.
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
    CSV.write(joinpath(surv_dir, "surviving.csv"), table)
    return surviving
end

# Agglomerative sub-clustering
# -----

# Full `N×N` pairwise CA-RMSD matrix for one set of conformers.
#
# `residues_list` is a vector of equal-length CA-only `PDBResidue` vectors (the members of one
# KMeans cluster, restricted to the ensemble-common set). Each pair is optimally superimposed with
# the package's MIToS superposition (the single source of truth for structural distance). The pair
# is re-superimposed rather than reusing the common-frame coordinates because the conformers are
# only aligned to the first conformer, so a common-frame RMSD would be an upper bound on the
# optimal pairwise RMSD. Returns a symmetric matrix with a zero diagonal.
function _pairwise_rmsd_matrix(residues_list::AbstractVector)
    n = length(residues_list)
    dist = zeros(Float64, n, n)
    for i = 1:n, j = (i+1):n
        _, _, rmsd = MIToS.PDB.superimpose(residues_list[i], residues_list[j])
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
# Columns: `Name`, `KMeans_K<k>` (the KMeans label, matching the kmeans table),
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

# Dimensionality reduction (MultivariateStats / Distances)
# -----

# Principal-component projection of a `D×N` feature matrix onto at most `ndims` axes.
#
# Columns are observations. Wraps `MultivariateStats.PCA` (`pratio=1.0`), returning the score
# matrix `predict` produces - `outdim×N`, where `outdim ≤ ndims` is the number of principal
# directions the data actually spans (rank-deficient inputs yield fewer rows; plotting pads to two
# dimensions via `_ensure_2d`). Deterministic given the input.
function _pca(features::AbstractMatrix, ndims::Integer = 2)
    model = MultivariateStats.fit(
        MultivariateStats.PCA,
        features;
        pratio = 1.0,
        maxoutdim = ndims,
    )
    return MultivariateStats.predict(model, features)
end

# Classical multidimensional scaling (PCoA) of an `N×N` distance matrix into at most `ndims` axes.
#
# Wraps `MultivariateStats.MDS` (classical MDS on a precomputed distance matrix), returning the
# `predict` coordinate matrix (one column per point) that recovers the original geometry up to
# rotation/reflection. An empty distance matrix yields an `ndims×0` layout. Deterministic given the
# input.
function _classical_mds(dist::AbstractMatrix, ndims::Integer = 2)
    n = size(dist, 1)
    n == 0 && return zeros(Float64, ndims, 0)
    model = MultivariateStats.fit(
        MultivariateStats.MDS,
        Matrix{Float64}(dist);
        distances = true,
        maxoutdim = ndims,
    )
    return MultivariateStats.predict(model)
end

# Full `N×N` euclidean distance matrix between the columns of a `D×N` feature matrix.
#
# Used to feed classical MDS for the embedding views; symmetric with a zero diagonal.
function _feature_distance_matrix(features::AbstractMatrix)
    return Distances.pairwise(Distances.Euclidean(), features; dims = 2)
end
# Centroid-to-centroid CA-RMSD between every occupied KMeans cluster.
#
# `features` is the `(3L)×N` aligned-coordinate matrix and `labels` the cluster labels.
# Returns `(cluster_ids, inter_rmsd)` where `cluster_ids` are the sorted occupied labels and
# `inter_rmsd` the matching `C×C` matrix of centroid euclidean distance divided by `sqrt(n_atoms)`
# (the CA-RMSD between cluster centroids), the query-path metric.
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
# tolerates a query whose residue numbering differs from the ensemble). `conformer_residues` are
# the already-read CA-only residues, in `names` order. Returns the cluster id. Raises a clear error
# when a name is not found and the value is not an existing file.
function _resolve_query_cluster(
    query::AbstractString,
    names::AbstractVector,
    labels::AbstractVector,
    conformer_residues::AbstractVector,
)
    pos = findfirst(==(query), names)
    pos !== nothing && return labels[pos]
    isfile(query) ||
        error("query '$query' is neither a conformer name nor an existing structure file")

    query_residues = _read_ca_residues(query)
    rmsds = _reference_rmsds(query_residues, conformer_residues)
    cluster, _mean = _assign_reference_cluster(rmsds, labels)
    return cluster
end

# Public entry point
# -----

"""
    cluster_conformers(system, ref1="", ref2="";
        ref_labels=("ref1", "ref2"), refs_dir="refs",
        data_root="/data/alphaconformers/pdb", out_dir="results",
        kmeans_k=30, subcluster_threshold=1.0, linkage=:average,
        score_table=nothing, score_col="cb-lddt", surviving_rel=25, surviving_frac=75,
        seed=42, make_plots=true, query=nothing, subdir=nothing)

Cluster an AlphaConformer prediction ensemble by shape, write the base KMeans tables, and
(when references are supplied) report each reference against the cluster it falls into.

It discovers the system's conformers, keeps the CA positions shared by all of them,
superimposes the coordinates, runs a seeded KMeans at `kmeans_k`, and writes the clustering
table and the intra-cluster RMSD table at the per-system output root. Up to two references may be
supplied as `ref1`/`ref2`; either may be given on its own. Each carries a user-chosen label from
`ref_labels` (default `"ref1"`/`"ref2"`). For each supplied reference it computes the
sequence-aware RMSD from every conformer, assigns the reference to the cluster whose conformers
are closest to it on average, and writes the per-reference RMSD table (a `RMSD_<label>` column per
reference) and the reference cluster-assignment table. When a DeepAccNet score table is supplied
it runs the optional score filter: for the configured relative/fraction cut it finds the surviving
clusters (a cluster survives when at most the fraction of its members fall below the relative score
cutoff) and writes the surviving-cluster table under a `deepaccnet/` folder; with no score table
this stage is skipped. It then always sub-clusters: within each KMeans cluster it builds the
pairwise RMSD matrix of the members, cuts an agglomerative tree at `subcluster_threshold`, and
writes that cluster's `members.csv` into its own `kmeans<C>/` folder, while a flat top-level
`agglomerative_assignments.csv` keeps the single-table view. When the score filter ran, only the
surviving KMeans clusters get a `kmeans<C>/` folder. Unless `make_plots` is disabled, it also
writes the figures for each level, adapting to the reference count (0/1/2). The query-path view is
added whenever a reference is given (it tracks the first reference) or an explicit `query` is
supplied.

# Arguments

- `system::String`: System identifier; names the input subdirectory under `data_root` and the
  output subdirectory under `out_dir`. Must not be blank after whitespace stripping.
- `ref1::String`: (Optional) first reference stem, e.g. `"1ABC_A"`. An empty/blank stem selects
  reference-free for this slot. Defaults to `""`.
- `ref2::String`: (Optional) second reference stem. May be given with or without `ref1`. Defaults
  to `""`.

# Keywords

- `ref_labels`: Two-element collection naming the references (used in the output columns and figure
  axes). Defaults to `("ref1", "ref2")`; set arbitrary strings to label the references, e.g.
  `("apo", "holo")`, `("monomer", "hexamer")` or `("-NaCl", "+NaCl")`.
- `refs_dir::String`: Directory holding the reference files. Each stem is resolved against this
  directory (`.pdb`/`.cif` and their `.gz` variants). Defaults to `"refs"`.
- `data_root::String`: Root directory holding one AlphaConformers output tree per `<system>`
  subfolder; predictions are read from the `models/` folders inside it. Defaults to
  `"/data/alphaconformers/pdb"`, matching where the package writes prediction outputs.
- `out_dir::String`: Root directory for output. Defaults to `"results"`.
- `kmeans_k::Int`: Number of clusters. Must be positive and at most the conformer count.
  Defaults to `30`.
- `subcluster_threshold::Float64`: Distance threshold (Å) at which the agglomerative tree is cut
  into sub-clusters. Must be positive. Defaults to `1.0`.
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
  falls back to the first supplied reference if any reference was given. So a run with references
  always writes the query-path figure tracking that reference, and passing `query` explicitly
  overrides this. In reference-free mode with no `query` the query-path figure is skipped.
- `subdir`: Optional relative subfolder under `data_root/<system>` to read conformers from (it may
  contain glob wildcards, e.g. `"cluster_*/af/predictions/sequences/models"`). When `nothing` (the
  default) the whole system directory is walked and only `models/` folders are read.

# Returns

A named tuple `(; clustering, cluster_rmsd, system_dir, reference_rmsd, reference_clusters,
deepaccnet_dir, surviving, hierarchical, cluster_dirs, figures)`:

- `clustering`: `DataFrame` with columns `Name` and `KMeans_K<k>` (labels in `1:k`).
- `cluster_rmsd`: `DataFrame` with columns `Method`, `Cluster`, `Size`, `Mean_RMSD_Angstrom`.
- `system_dir`: path to the per-system output root holding the flat tables, the system-level
  figures, the `kmeans<C>/` folders and the optional `deepaccnet/` folder.
- `reference_rmsd`: `nothing` in reference-free mode; otherwise a `DataFrame` with a `Name`
  column and a `RMSD_<label>` column per supplied reference, one row per conformer.
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

- `ErrorException`: if `system` is blank, if `kmeans_k <= 0` or `subcluster_threshold <= 0`, if a
  reference stem cannot be resolved under `refs_dir`, if no conformers are found, if the conformers
  share no common CA position, if `kmeans_k` exceeds the conformer count, or, when a score table is
  supplied, if it is missing the `Name`/`score_col` column, if it matches no conformer, if
  `surviving_rel` or `surviving_frac` is outside `0..100`, or if a `query` is given that is neither
  a known conformer name nor an existing structure file.

# Behavior

1. Validate inputs and count the references supplied (0, 1 or 2; either positional may be given on
   its own).
2. Discover the prediction structures under `data_root/<system>` - by default only the files in
   `models/` folders (the AlphaConformers prediction location), so input templates are skipped, or
   the explicit `subdir` when given - read their CA traces, and keep the CA positions shared by all
   conformers via sequence alignment.
3. Superimpose every conformer onto the first and run seeded KMeans at `kmeans_k`.
4. Compute the mean intra-cluster RMSD to each cluster centroid and write
   `aligned_clustering_results.csv` and `aligned_cluster_rmsd.csv` at `out_dir/<system>/`.
5. With one or two references, resolve each stem against `refs_dir`, compute the sequence-aware
   RMSD from every conformer to each reference, assign each reference to the cluster with the
   lowest mean RMSD, and write `dist_external_rmsds.csv` and `reference_clusters.csv` at
   `out_dir/<system>/`.
6. When a score table is supplied, align its scores to the conformers (dropping unscored ones),
   find the surviving clusters for the configured `surviving_rel`/`surviving_frac` cut, and write
   `surviving.csv` under `out_dir/<system>/deepaccnet/surviving/`. Skipped
   when no score table is given.
7. For each occupied (or, when the filter ran, surviving) KMeans cluster `C`, build the pairwise
   RMSD matrix of its members, cut the agglomerative tree at `subcluster_threshold` with the chosen
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
    ref1::String = "",
    ref2::String = "";
    ref_labels = ("ref1", "ref2"),
    refs_dir::String = "refs",
    data_root::String = "/data/alphaconformers/pdb",
    out_dir::String = "results",
    kmeans_k::Int = 30,
    subcluster_threshold::Float64 = 1.0,
    linkage::Symbol = :average,
    score_table = nothing,
    score_col::String = "cb-lddt",
    surviving_rel::Int = 25,
    surviving_frac::Int = 75,
    seed::Int = 42,
    make_plots::Bool = true,
    query::Union{Nothing,AbstractString} = nothing,
    subdir::Union{Nothing,AbstractString} = nothing,
)
    # Validate and normalize inputs.
    stripped_system = strip(system)
    isempty(stripped_system) && error("system must not be blank")

    kmeans_k > 0 || error("kmeans_k must be positive")
    subcluster_threshold > 0 || error("subcluster_threshold must be positive")

    ref1_normalized = _normalize_reference(ref1)
    ref2_normalized = _normalize_reference(ref2)

    # Number of references supplied (0, 1 or 2); either positional may be given on its own.
    ref_mode = _reference_mode(ref1_normalized, ref2_normalized)

    # Discover the ensemble's conformers and read their CA-only residues.
    conformers = _discover_conformers(data_root, String(stripped_system); subdir = subdir)
    names = first.(conformers)
    conformer_residues = [_read_ca_residues(path) for (_name, path) in conformers]

    if kmeans_k > length(conformers)
        error(
            "kmeans_k ($kmeans_k) exceeds the number of conformers " *
            "($(length(conformers))); choose a smaller kmeans_k.",
        )
    end

    # Build the ensemble-wide common-Cα set by sequence alignment, then superimpose.
    common_residues = _ensemble_common_residues(conformer_residues)
    aligned = _superimpose_common(common_residues)
    features = _feature_matrix(aligned)
    n_atoms = length(first(common_residues))

    # Seeded KMeans and intra-cluster RMSD.
    labels = _kmeans_cluster(features, kmeans_k; seed = seed)
    rmsd = _intra_cluster_rmsd(features, labels, n_atoms)

    # Per-system output root: the top-level flat tables and the system-level figures live here,
    # alongside the per-cluster `kmeans<C>/` folders and the optional `deepaccnet/` folder.
    system_dir = joinpath(out_dir, String(stripped_system))
    clustering, cluster_rmsd =
        _write_kmeans_tables(system_dir, names, labels, kmeans_k, rmsd)

    # Reference handling: sequence-aware RMSD from each supplied reference to every conformer and
    # the cluster it falls into (lowest mean RMSD). Each reference carries its user-chosen label
    # (default `ref1`/`ref2`). Skipped entirely in reference-free mode.
    reference_rmsd = nothing
    reference_clusters = nothing
    rmsd_ref1 = Float64[]
    rmsd_ref2 = Float64[]
    present_labels = String[]
    if ref_mode > 0
        refs = Tuple{String,String,Vector{Float64}}[]
        for (label, stem) in zip(ref_labels, (ref1_normalized, ref2_normalized))
            stem === nothing && continue
            rmsds = _reference_rmsds(
                _read_ca_residues(_resolve_reference(refs_dir, stem)),
                conformer_residues,
            )
            push!(refs, (String(label), stem, rmsds))
        end

        reference_rmsd, reference_clusters =
            _write_reference_tables(system_dir, names, labels, refs)
        present_labels = [r[1] for r in refs]
        rmsd_ref1 = refs[1][3]
        ref_mode == 2 && (rmsd_ref2 = refs[2][3])
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
    # RMSD matrix of its members, cut the agglomerative tree at `subcluster_threshold`, and record a
    # per-conformer sub-cluster label. When the score filter ran, only the scored conformers in
    # the configured surviving clusters are sub-clustered; otherwise the full KMeans clustering.
    # Each cluster's outputs land in its own `kmeans<C>/` folder as a flat `members.csv` (no leaf
    # sub-cluster folders); a flat top-level `agglomerative_assignments.csv` keeps the single-table
    # view.
    if score_table === nothing
        work_idx = collect(eachindex(names))
    else
        keep = Set(surviving)
        work_idx = [i for i in scored_idx if labels[i] in keep]
    end
    work_names = names[work_idx]
    work_labels = labels[work_idx]
    work_residues = common_residues[work_idx]

    sub_labels = zeros(Int, length(work_idx))
    cluster_results = NamedTuple[]
    cluster_dirs = String[]
    for lbl in sort(unique(work_labels))
        idx = findall(==(lbl), work_labels)
        dist = _pairwise_rmsd_matrix(work_residues[idx])
        subs = _agglomerative_subcluster(dist, subcluster_threshold, linkage)
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

    # Query-path inputs (C4). The query defaults to the first supplied reference, so the query-path
    # view tracks that reference without an extra argument. An explicit `query` keyword overrides
    # this default. Centroid distances are computed over the full KMeans clustering so the ranking
    # covers every occupied cluster. For the defaulted reference query the cluster is that
    # reference's already-computed assigned cluster (row 1 of `reference_clusters`), avoiding a
    # second RMSD pass.
    query_cluster = nothing
    cluster_ids = Int[]
    inter_rmsd = zeros(0, 0)
    if query !== nothing
        cluster_ids, inter_rmsd = _cluster_centroid_rmsd(features, labels, n_atoms)
        query_cluster = _resolve_query_cluster(query, names, labels, conformer_residues)
    elseif ref_mode > 0
        cluster_ids, inter_rmsd = _cluster_centroid_rmsd(features, labels, n_atoms)
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
            rmsd_ref1,
            rmsd_ref2,
            ref_labels = present_labels,
            reference_clusters,
            scored_idx,
            scores,
            surviving,
            cluster_results,
            linkage,
            threshold = subcluster_threshold,
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
