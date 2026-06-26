using TestItems

@testitem "Conformer clustering input validation" begin
    # Test _normalize_reference helper.
    @test AlphaConformers._normalize_reference("") === nothing
    @test AlphaConformers._normalize_reference("   ") === nothing
    @test AlphaConformers._normalize_reference("1ABC_A") == "1ABC_A"

    # Test _reference_mode helper: just the count of references supplied; either slot may be given
    # on its own (no required ordering).
    @test AlphaConformers._reference_mode(nothing, nothing) == 0  # 0-ref
    @test AlphaConformers._reference_mode("1ABC_A", nothing) == 1  # 1-ref (first slot)
    @test AlphaConformers._reference_mode(nothing, "1DEF_B") == 1  # 1-ref (second slot only)
    @test AlphaConformers._reference_mode("1ABC_A", "1DEF_B") == 2  # 2-ref

    # Test cluster_conformers validation errors.
    @test_throws ErrorException cluster_conformers("")  # blank system
    @test_throws ErrorException cluster_conformers("   ")  # whitespace-only system
    @test_throws ErrorException cluster_conformers("sys1"; kmeans_k = 0)  # invalid k
    @test_throws ErrorException cluster_conformers("sys1"; kmeans_k = -5)  # negative k
    @test_throws ErrorException cluster_conformers("sys1"; subcluster_threshold = 0.0)
    @test_throws ErrorException cluster_conformers("sys1"; subcluster_threshold = -1.5)
    # A second-slot reference on its own is allowed (no asymmetry); this still errors only because
    # there is no discoverable data at the default data_root.
    @test_throws ErrorException cluster_conformers("sys1", "", "1ABC_A")

    # Reference modes run now; with no discoverable conformers the run still errors clearly.
    @test_throws ErrorException cluster_conformers(
        "sys1",
        "1ABC_A";
        data_root = mktempdir(),
    )
    @test_throws ErrorException cluster_conformers(
        "sys1",
        "1ABC_A",
        "1DEF_B";
        data_root = mktempdir(),
    )

    # A reference-free call with no discoverable data errors clearly.
    @test_throws ErrorException cluster_conformers("sys1"; data_root = mktempdir())
end

@testitem "Conformer clustering KMeans walking skeleton" begin
    import CSV
    import DataFrames
    import Random

    # Build a tiny, deterministic conformer ensemble in a temp data root: three shape
    # groups (A/B/C) of three replicates each, ten CA residues, differing by a
    # rotation-invariant progressive bend so the groups stay distinct after superposition.
    function _write_fixture(data_root, system)
        models = joinpath(data_root, system, "models")
        mkpath(models)
        residues = 10
        base = [(3.8 * (i - 1), 0.0, 0.0) for i = 1:residues]
        for (group, bend) in ["A" => 0.0, "B" => 5.0, "C" => 12.0]
            for rep = 1:3
                rng = Random.MersenneTwister(abs(hash((group, rep))) % 100_000)
                path = joinpath(models, "conf_$(group)$(rep).pdb")
                open(path, "w") do io
                    for (i, b) in enumerate(base)
                        yb =
                            i > residues ÷ 2 ? bend * (i - residues ÷ 2) / (residues ÷ 2) :
                            0.0
                        jitter = 0.05 .* (Random.rand(rng, 3) .- 0.5)
                        x, y, z = b[1] + jitter[1], b[2] + yb + jitter[2], b[3] + jitter[3]
                        line = string(
                            "ATOM  ",
                            lpad(i, 5),
                            "  CA  ALA A",
                            lpad(i, 4),
                            "    ",
                            lpad(string(round(x, digits = 3)), 8),
                            lpad(string(round(y, digits = 3)), 8),
                            lpad(string(round(z, digits = 3)), 8),
                            "  1.00  0.00           C",
                        )
                        println(io, rpad(line, 80))
                    end
                    println(io, "END")
                end
            end
        end
        return data_root
    end

    system = "toy_system"
    k = 3
    n_conformers = 9
    fixtures = _write_fixture(mktempdir(), system)

    out = mktempdir()
    result = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = out,
        kmeans_k = k,
        seed = 42,
    )

    # The per-system output root and both top-level tables are written.
    @test result.system_dir == joinpath(out, system)
    @test isdir(result.system_dir)
    clustering_csv = joinpath(result.system_dir, "aligned_clustering_results.csv")
    rmsd_csv = joinpath(result.system_dir, "aligned_cluster_rmsd.csv")
    @test isfile(clustering_csv)
    @test isfile(rmsd_csv)

    # Clustering table schema and invariants (read back from disk).
    clustering = CSV.read(clustering_csv, DataFrames.DataFrame)
    label_col = "KMeans_K$k"
    @test names(clustering) == ["Name", label_col]
    @test DataFrames.nrow(clustering) == n_conformers          # every conformer is labelled
    labels = clustering[!, label_col]
    @test !any(ismissing, labels)
    @test all(1 .<= labels .<= k)                              # labels fall in 1:k
    @test length(unique(clustering.Name)) == n_conformers      # one row per conformer

    # Intra-cluster RMSD table schema and invariants.
    cluster_rmsd = CSV.read(rmsd_csv, DataFrames.DataFrame)
    @test names(cluster_rmsd) == ["Method", "Cluster", "Size", "Mean_RMSD_Angstrom"]
    @test all(cluster_rmsd.Method .== label_col)
    @test all(1 .<= cluster_rmsd.Cluster .<= k)
    @test sum(cluster_rmsd.Size) == n_conformers               # sizes cover all conformers
    @test all(cluster_rmsd.Mean_RMSD_Angstrom .>= 0)

    # The returned DataFrames match what was written.
    @test result.clustering == clustering
    @test result.cluster_rmsd == cluster_rmsd

    # Reproducibility: same inputs and seed yield identical labels.
    again = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = mktempdir(),
        kmeans_k = k,
        seed = 42,
    )
    @test again.clustering[!, label_col] == labels

    # kmeans_k larger than the conformer count is a clear error.
    @test_throws ErrorException cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = mktempdir(),
        kmeans_k = n_conformers + 1,
    )
end

@testitem "Conformer clustering helper seams" begin
    import MIToS

    # Build a vector of CA-only poly-ALA PDBResidues from an L×3 coordinate matrix (the test
    # fixtures' counterpart to a real structure read).
    function make_ca_residues(coords)
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

    # _intra_cluster_rmsd: centroid metric is exact; singletons get 0.0.
    # Columns are observations; here one CA atom (3 coordinates), so n_atoms = 1.
    spread = Float64[0 2; 0 0; 0 0]                 # two members at x=0 and x=2
    @test AlphaConformers._intra_cluster_rmsd(spread, [7, 7], 1)[7] ≈ 1.0  # |±1| → mean 1.0

    mixed = Float64[0 0 5; 0 0 0; 0 0 0]            # cluster 1 = two identical; cluster 2 singleton
    rmsd = AlphaConformers._intra_cluster_rmsd(mixed, [1, 1, 2], 1)
    @test rmsd[1] == 0.0                             # identical members
    @test rmsd[2] == 0.0                             # singleton cluster branch

    # _read_ca_residues reads CA traces from .cif as well as .pdb, and CAmatrix recovers them.
    tmp = mktempdir()
    residues = make_ca_residues(Float64[0 0 0; 3.8 0 0; 7.6 0 0])
    cif = joinpath(tmp, "ca_only.cif")
    MIToS.PDB.write_file(cif, residues, MIToS.PDB.MMCIFFile)
    got = AlphaConformers._read_ca_residues(cif)
    coords = MIToS.PDB.CAmatrix(got)
    @test length(got) == 3
    @test size(coords) == (3, 3)
    @test coords[2, 1] ≈ 3.8
end

@testitem "Conformer clustering reference resolution" begin
    import MIToS

    refs = mktempdir()
    touch(joinpath(refs, "1ABC_A.pdb"))    # exact .pdb match
    touch(joinpath(refs, "9XYZ_B.cif"))    # only .cif present
    touch(joinpath(refs, "2def_C.pdb"))    # lowercase-id file for stem 2DEF_C
    touch(joinpath(refs, "3GZP_A.cif.gz")) # only a gzip-compressed reference present

    # Exact stem, .pdb preferred.
    @test AlphaConformers._resolve_reference(refs, "1ABC_A") == joinpath(refs, "1ABC_A.pdb")
    # `.cif` fallback when no `.pdb` exists.
    @test AlphaConformers._resolve_reference(refs, "9XYZ_B") == joinpath(refs, "9XYZ_B.cif")
    # A compressed `.cif.gz` reference resolves (MIToS decompresses `.gz` on read).
    @test AlphaConformers._resolve_reference(refs, "3GZP_A") ==
          joinpath(refs, "3GZP_A.cif.gz")
    # `_structure_format` looks past a trailing `.gz` to pick the MIToS format.
    @test AlphaConformers._structure_format("x.pdb") == MIToS.PDB.PDBFile
    @test AlphaConformers._structure_format("x.cif") == MIToS.PDB.MMCIFFile
    @test AlphaConformers._structure_format("x.pdb.gz") == MIToS.PDB.PDBFile
    @test AlphaConformers._structure_format("x.cif.gz") == MIToS.PDB.MMCIFFile
    # Case-insensitive lowercase-id fallback that keeps the chain letter. Compared by file
    # identity rather than path string, so a case-insensitive filesystem (which resolves the
    # exact-case stem to the same on-disk file) passes too.
    @test Base.Filesystem.samefile(
        AlphaConformers._resolve_reference(refs, "2DEF_C"),
        joinpath(refs, "2def_C.pdb"),
    )
    # Unresolvable stem raises a clear error.
    @test_throws ErrorException AlphaConformers._resolve_reference(refs, "0NON_E")
end

@testitem "Conformer clustering reference RMSD seams" begin
    import MIToS

    function make_ca_residues(coords)
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

    # CA-only poly-ALA structures built from coordinate matrices (a small 3-D shape so the
    # superposition is non-degenerate).
    coords = Float64[0 0 0; 3.8 0 0; 3.8 3.8 0; 0 3.8 0; 1.9 1.9 3.0]
    ref = make_ca_residues(coords)

    # Self-RMSD ≈ 0.
    self = make_ca_residues(coords)
    @test AlphaConformers._reference_rmsds(ref, [self])[1] ≈ 0.0 atol = 1e-6

    # Pure translation is removed by the superposition → ≈ 0.
    translated = make_ca_residues(coords .+ [5.0 -3.0 2.0])
    @test AlphaConformers._reference_rmsds(ref, [translated])[1] ≈ 0.0 atol = 1e-6

    # Known pair: scaling about the centroid by `c` cannot be undone by a rigid move, so the
    # RMSD equals |1 - c| times the radius of gyration (optimal rotation is the identity).
    c = 1.1
    natoms = size(coords, 1)
    centroid = sum(coords; dims = 1) ./ natoms
    scaled = centroid .+ c .* (coords .- centroid)
    conf = make_ca_residues(scaled)
    rgyr = sqrt(sum(abs2, coords .- centroid) / natoms)
    @test AlphaConformers._reference_rmsds(ref, [conf])[1] ≈ abs(1 - c) * rgyr atol = 1e-4

    # A failed alignment (no shared residues / nothing) surfaces as NaN, not an error.
    @test isnan(AlphaConformers._reference_rmsds(ref, [MIToS.PDB.PDBResidue[]])[1])

    # _read_ca_residues reads a `.cif` reference (CA-only) as well as `.pdb`.
    cif = joinpath(mktempdir(), "ref.cif")
    MIToS.PDB.write_file(cif, ref, MIToS.PDB.MMCIFFile)
    cif_res = AlphaConformers._read_ca_residues(cif)
    @test length(cif_res) == size(coords, 1)
    @test all(r -> !isempty(MIToS.PDB.findatoms(r, "CA")), cif_res)

    # _assign_reference_cluster: lowest mean RMSD wins; NaN ignored; all-NaN cluster skipped.
    rmsds = [0.1, 0.2, 5.0, 5.0, NaN, 0.05]
    labels = [1, 1, 2, 2, 3, 3]
    cluster, mean_rmsd = AlphaConformers._assign_reference_cluster(rmsds, labels)
    @test cluster == 3                 # cluster 3's only finite value (0.05) is the smallest mean
    @test mean_rmsd ≈ 0.05
    @test_throws ErrorException AlphaConformers._assign_reference_cluster(
        [NaN, NaN],
        [1, 2],
    )

    # Modified residues are handled by MIToS' modelled_sequences (used inside structural_alignment),
    # so the pipeline needs no bespoke modified-residue table: a selenomethionine (MSE) maps to the
    # standard methionine letter in the extracted sequence.
    modres = make_ca_residues(Float64[0 0 0; 3.8 0 0; 7.6 0 0])
    old = modres[2].id
    modres[2].id = MIToS.PDB.PDBResidueIdentifier(
        old.PDBe_number,
        old.number,
        "MSE",
        "ATOM",
        old.model,
        old.chain,
    )
    seq = string(first(values(MIToS.PDB.modelled_sequences(modres))))
    @test seq[2] == 'M'
end

@testitem "Conformer clustering with references" begin
    import CSV
    import DataFrames
    import Random

    # Same deterministic ensemble as the walking-skeleton test (three shape groups).
    function _write_fixture(data_root, system)
        models = joinpath(data_root, system, "models")
        mkpath(models)
        residues = 10
        base = [(3.8 * (i - 1), 0.0, 0.0) for i = 1:residues]
        for (group, bend) in ["A" => 0.0, "B" => 5.0, "C" => 12.0]
            for rep = 1:3
                rng = Random.MersenneTwister(abs(hash((group, rep))) % 100_000)
                path = joinpath(models, "conf_$(group)$(rep).pdb")
                open(path, "w") do io
                    for (i, b) in enumerate(base)
                        yb =
                            i > residues ÷ 2 ? bend * (i - residues ÷ 2) / (residues ÷ 2) :
                            0.0
                        jitter = 0.05 .* (Random.rand(rng, 3) .- 0.5)
                        x, y, z = b[1] + jitter[1], b[2] + yb + jitter[2], b[3] + jitter[3]
                        line = string(
                            "ATOM  ",
                            lpad(i, 5),
                            "  CA  ALA A",
                            lpad(i, 4),
                            "    ",
                            lpad(string(round(x, digits = 3)), 8),
                            lpad(string(round(y, digits = 3)), 8),
                            lpad(string(round(z, digits = 3)), 8),
                            "  1.00  0.00           C",
                        )
                        println(io, rpad(line, 80))
                    end
                    println(io, "END")
                end
            end
        end
        return data_root
    end

    system = "toy_system"
    k = 3
    n_conformers = 9
    fixtures = _write_fixture(mktempdir(), system)
    models = joinpath(fixtures, system, "models")

    # References are pre-trimmed single-chain files; reuse two conformers as the apo/holo
    # references so each reference is identical to one ensemble member.
    refs = mktempdir()
    cp(joinpath(models, "conf_A1.pdb"), joinpath(refs, "REFAP_A.pdb"))
    cp(joinpath(models, "conf_C1.pdb"), joinpath(refs, "REFHO_B.pdb"))

    # Baseline (reference-free) labels, to know which cluster conf_A1 / conf_C1 land in.
    base = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = mktempdir(),
        kmeans_k = k,
        seed = 42,
    )
    label_col = "KMeans_K$k"
    label_of(name) = base.clustering[findfirst(==(name), base.clustering.Name), label_col]

    # 1-reference mode (apo only).
    out1 = mktempdir()
    one = cluster_conformers(
        system,
        "REFAP_A";
        refs_dir = refs,
        data_root = fixtures,
        out_dir = out1,
        kmeans_k = k,
        seed = 42,
    )

    rmsd_csv = joinpath(out1, system, "dist_external_rmsds.csv")
    refclu_csv = joinpath(one.system_dir, "reference_clusters.csv")
    @test isfile(rmsd_csv)
    @test isfile(refclu_csv)

    rmsd_tbl = CSV.read(rmsd_csv, DataFrames.DataFrame)
    # Default labels are the generic ref1/ref2, so the column follows the first label.
    @test names(rmsd_tbl) == ["Name", "RMSD_ref1"]         # one column for the single reference
    @test DataFrames.nrow(rmsd_tbl) == n_conformers        # one row per conformer
    @test one.reference_rmsd == rmsd_tbl

    refclu = CSV.read(refclu_csv, DataFrames.DataFrame)
    @test names(refclu) == ["Reference", "Stem", "Cluster", "Mean_RMSD_Angstrom"]
    @test refclu.Reference == ["ref1"]
    @test refclu.Stem == ["REFAP_A"]
    @test all(1 .<= refclu.Cluster .<= k)                 # a valid cluster id
    @test one.reference_clusters == refclu
    # The first reference is conf_A1, so it is assigned to conf_A1's cluster.
    @test refclu.Cluster[1] == label_of("models/conf_A1.pdb")
    # That conformer's own RMSD to the reference is ~0.
    ref1_row = rmsd_tbl[findfirst(==("models/conf_A1.pdb"), rmsd_tbl.Name), :]
    @test ref1_row.RMSD_ref1 ≈ 0.0 atol = 1e-4

    # Custom `ref_labels` flow through to the column and assignment labels.
    out_lbl = mktempdir()
    labelled = cluster_conformers(
        system,
        "REFAP_A",
        "REFHO_B";
        ref_labels = ("apo", "holo"),
        refs_dir = refs,
        data_root = fixtures,
        out_dir = out_lbl,
        kmeans_k = k,
        seed = 42,
    )
    @test names(labelled.reference_rmsd) == ["Name", "RMSD_apo", "RMSD_holo"]
    @test labelled.reference_clusters.Reference == ["apo", "holo"]

    # 0-reference mode leaves the reference fields empty.
    @test base.reference_rmsd === nothing
    @test base.reference_clusters === nothing

    # 2-reference mode (default ref1/ref2 labels).
    out2 = mktempdir()
    two = cluster_conformers(
        system,
        "REFAP_A",
        "REFHO_B";
        refs_dir = refs,
        data_root = fixtures,
        out_dir = out2,
        kmeans_k = k,
        seed = 42,
    )
    rmsd2 =
        CSV.read(joinpath(out2, system, "dist_external_rmsds.csv"), DataFrames.DataFrame)
    @test names(rmsd2) == ["Name", "RMSD_ref1", "RMSD_ref2"]
    refclu2 = two.reference_clusters
    @test refclu2.Reference == ["ref1", "ref2"]
    @test all(1 .<= refclu2.Cluster .<= k)
    # apo == conf_A1, holo == conf_C1, each assigned to its own conformer's cluster.
    @test refclu2.Cluster[1] == label_of("models/conf_A1.pdb")
    @test refclu2.Cluster[2] == label_of("models/conf_C1.pdb")

    # Unresolvable reference stem raises a clear error.
    @test_throws ErrorException cluster_conformers(
        system,
        "NOPE_A";
        refs_dir = refs,
        data_root = fixtures,
        out_dir = mktempdir(),
        kmeans_k = k,
    )
end

@testitem "Conformer clustering agglomerative cut" begin
    # Relabel a label vector by first occurrence so the assertion is exact regardless of how
    # cutree numbers the sub-clusters internally.
    function canonical(labels)
        seen = Dict{Int,Int}()
        out = similar(labels)
        for (i, l) in enumerate(labels)
            out[i] = get!(seen, l, length(seen) + 1)
        end
        return out
    end

    # Two tight pairs ({1,2} and {3,4}) sitting far apart.
    dist = [
        0.0 0.5 5.0 5.0
        0.5 0.0 5.0 5.0
        5.0 5.0 0.0 0.5
        5.0 5.0 0.5 0.0
    ]

    # Cut at 1.0 Å: the two close pairs split from each other.
    @test canonical(AlphaConformers._agglomerative_subcluster(dist, 1.0, :average)) ==
          [1, 1, 2, 2]
    # A high threshold collapses everything into one sub-cluster.
    @test AlphaConformers._agglomerative_subcluster(dist, 10.0, :average) == [1, 1, 1, 1]
    # A tiny threshold makes every member its own sub-cluster.
    @test canonical(AlphaConformers._agglomerative_subcluster(dist, 0.1, :average)) ==
          [1, 2, 3, 4]

    # Linkage is configurable; :single and :complete both give the same two-pair split here.
    @test canonical(AlphaConformers._agglomerative_subcluster(dist, 1.0, :single)) ==
          [1, 1, 2, 2]
    @test canonical(AlphaConformers._agglomerative_subcluster(dist, 1.0, :complete)) ==
          [1, 1, 2, 2]

    # Degenerate cluster sizes: a single member is its own sub-cluster; empty yields nothing.
    @test AlphaConformers._agglomerative_subcluster(zeros(1, 1), 1.0, :average) == [1]
    @test AlphaConformers._agglomerative_subcluster(zeros(0, 0), 1.0, :average) == Int[]
end

@testitem "Conformer clustering pairwise RMSD matrix" begin
    import MIToS

    function make_ca_residues(coords)
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

    # Two identical squares and a translated copy (as CA-only residue vectors): the matrix is
    # symmetric, zero on the diagonal, and a pure translation is removed by the superposition
    # (off-diagonal ≈ 0).
    square = make_ca_residues(Float64[0 0 0; 3.8 0 0; 3.8 3.8 0; 0 3.8 0])
    translated =
        make_ca_residues(Float64[0 0 0; 3.8 0 0; 3.8 3.8 0; 0 3.8 0] .+ [5.0 -2.0 1.0])
    dist = AlphaConformers._pairwise_rmsd_matrix([square, deepcopy(square), translated])

    @test size(dist) == (3, 3)
    @test dist == dist'                       # symmetric
    @test all(dist[i, i] == 0.0 for i = 1:3)  # zero diagonal
    @test dist[1, 2] ≈ 0.0 atol = 1e-6        # identical members
    @test dist[1, 3] ≈ 0.0 atol = 1e-6        # translation removed by superposition
end

@testitem "Conformer clustering hierarchical assignments" begin
    import CSV
    import DataFrames
    import Random

    # Same deterministic ensemble as the walking-skeleton test (three shape groups).
    function _write_fixture(data_root, system)
        models = joinpath(data_root, system, "models")
        mkpath(models)
        residues = 10
        base = [(3.8 * (i - 1), 0.0, 0.0) for i = 1:residues]
        for (group, bend) in ["A" => 0.0, "B" => 5.0, "C" => 12.0]
            for rep = 1:3
                rng = Random.MersenneTwister(abs(hash((group, rep))) % 100_000)
                path = joinpath(models, "conf_$(group)$(rep).pdb")
                open(path, "w") do io
                    for (i, b) in enumerate(base)
                        yb =
                            i > residues ÷ 2 ? bend * (i - residues ÷ 2) / (residues ÷ 2) :
                            0.0
                        jitter = 0.05 .* (Random.rand(rng, 3) .- 0.5)
                        x, y, z = b[1] + jitter[1], b[2] + yb + jitter[2], b[3] + jitter[3]
                        line = string(
                            "ATOM  ",
                            lpad(i, 5),
                            "  CA  ALA A",
                            lpad(i, 4),
                            "    ",
                            lpad(string(round(x, digits = 3)), 8),
                            lpad(string(round(y, digits = 3)), 8),
                            lpad(string(round(z, digits = 3)), 8),
                            "  1.00  0.00           C",
                        )
                        println(io, rpad(line, 80))
                    end
                    println(io, "END")
                end
            end
        end
        return data_root
    end

    system = "toy_system"
    k = 3
    n_conformers = 9
    fixtures = _write_fixture(mktempdir(), system)

    out = mktempdir()
    result = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = out,
        kmeans_k = k,
        seed = 42,
    )

    # The flat top-level assignment table is written at the per-system root.
    assign_csv = joinpath(result.system_dir, "agglomerative_assignments.csv")
    @test isfile(assign_csv)

    assignments = CSV.read(assign_csv, DataFrames.DataFrame)
    label_col = "KMeans_K$k"
    @test names(assignments) == ["Name", label_col, "Sub_Cluster", "Mini_Cluster"]
    @test DataFrames.nrow(assignments) == n_conformers       # covers every conformer
    @test length(unique(assignments.Name)) == n_conformers   # one row per conformer

    # Each occupied KMeans cluster gets its own kmeans<C>/ folder holding a flat members.csv
    # (this cluster's members with their sub-cluster labels) and no leaf sub-folders.
    @test !isempty(result.cluster_dirs)
    for cluster_dir in result.cluster_dirs
        @test isdir(cluster_dir)
        members_csv = joinpath(cluster_dir, "members.csv")
        @test isfile(members_csv)
        members = CSV.read(members_csv, DataFrames.DataFrame)
        @test names(members) == ["Name", label_col, "Sub_Cluster", "Mini_Cluster"]
        @test length(unique(members[!, label_col])) == 1   # one KMeans cluster per folder
        # No nested per-sub-cluster folders remain.
        @test !any(isdir, joinpath.(cluster_dir, readdir(cluster_dir)))
    end
    # The occupied clusters each map to a folder, and members.csv files partition the ensemble.
    @test length(result.cluster_dirs) == length(unique(result.clustering[!, label_col]))
    @test sum(
        DataFrames.nrow(CSV.read(joinpath(d, "members.csv"), DataFrames.DataFrame)) for
        d in result.cluster_dirs
    ) == n_conformers

    # Every member of every cluster receives a sub-cluster label (≥ 1, no missing).
    @test !any(ismissing, assignments.Sub_Cluster)
    @test all(assignments.Sub_Cluster .>= 1)
    @test !any(ismissing, assignments.Mini_Cluster)

    # The KMeans labels match the base clustering, and the mini-cluster name is consistent.
    @test assignments[!, label_col] == result.clustering[!, label_col]
    expected_mini = [
        "c$(l)_s$(s)" for (l, s) in zip(assignments[!, label_col], assignments.Sub_Cluster)
    ]
    @test assignments.Mini_Cluster == expected_mini

    # Within each KMeans cluster, sub-cluster labels start at 1 and are contiguous.
    for lbl in unique(assignments[!, label_col])
        subs = assignments.Sub_Cluster[assignments[!, label_col] .== lbl]
        @test sort(unique(subs)) == collect(1:length(unique(subs)))
    end

    # The returned DataFrame matches what was written.
    @test result.hierarchical == assignments
end

@testitem "Conformer clustering surviving clusters" begin
    # Crafted scores + cutoffs yield an exact surviving set. Two clusters of three members each.
    # Scores: cluster 0 has two low and one high; cluster 1 has one low and two high.
    labels = [0, 0, 0, 1, 1, 1]
    scores = [0.1, 0.2, 0.9, 0.8, 0.85, 0.95]

    # The median (rel 50) cutoff is 0.825: cluster 0 has 2/3 of its members below it, cluster 1
    # only 1/3. So cluster 1 always survives, and cluster 0 survives only once the allowed
    # fraction reaches 75%.
    @test AlphaConformers._surviving_clusters(labels, scores, 50, 50) == [1]
    @test AlphaConformers._surviving_clusters(labels, scores, 50, 60) == [1]
    @test AlphaConformers._surviving_clusters(labels, scores, 50, 75) == [0, 1]

    # With nobody below the cutoff (rel 0 → the minimum, strict `<`), every cluster survives.
    @test AlphaConformers._surviving_clusters(labels, scores, 0, 50) == [0, 1]
end

@testitem "Conformer clustering score alignment" begin
    import DataFrames

    table = DataFrames.DataFrame("Name" => ["a", "b", "c"], "cb-lddt" => [0.5, 0.9, 0.1])

    # Names are matched in order; an unscored conformer (`x`) is dropped by the inner join.
    idx, scores = AlphaConformers._align_scores(table, ["a", "c", "x"])
    @test idx == [1, 2]
    @test scores == [0.5, 0.1]

    idx2, scores2 = AlphaConformers._align_scores(table, ["a", "b", "c", "x"])
    @test idx2 == [1, 2, 3]
    @test scores2 == [0.5, 0.9, 0.1]

    # A missing `Name` column, a missing score column, and no matching conformer all error.
    no_name = DataFrames.DataFrame("Id" => ["a"], "cb-lddt" => [0.1])
    @test_throws ErrorException AlphaConformers._align_scores(no_name, ["a"])
    no_score = DataFrames.DataFrame("Name" => ["a"], "other" => [0.1])
    @test_throws ErrorException AlphaConformers._align_scores(
        no_score,
        ["a"];
        score_col = "cb-lddt",
    )
    @test_throws ErrorException AlphaConformers._align_scores(table, ["zzz"])
end

@testitem "Conformer clustering DeepAccNet score filter" begin
    import CSV
    import DataFrames
    import Random

    function _write_fixture(data_root, system)
        models = joinpath(data_root, system, "models")
        mkpath(models)
        residues = 10
        base = [(3.8 * (i - 1), 0.0, 0.0) for i = 1:residues]
        for (group, bend) in ["A" => 0.0, "B" => 5.0, "C" => 12.0]
            for rep = 1:3
                rng = Random.MersenneTwister(abs(hash((group, rep))) % 100_000)
                path = joinpath(models, "conf_$(group)$(rep).pdb")
                open(path, "w") do io
                    for (i, b) in enumerate(base)
                        yb =
                            i > residues ÷ 2 ? bend * (i - residues ÷ 2) / (residues ÷ 2) :
                            0.0
                        jitter = 0.05 .* (Random.rand(rng, 3) .- 0.5)
                        x, y, z = b[1] + jitter[1], b[2] + yb + jitter[2], b[3] + jitter[3]
                        line = string(
                            "ATOM  ",
                            lpad(i, 5),
                            "  CA  ALA A",
                            lpad(i, 4),
                            "    ",
                            lpad(string(round(x, digits = 3)), 8),
                            lpad(string(round(y, digits = 3)), 8),
                            lpad(string(round(z, digits = 3)), 8),
                            "  1.00  0.00           C",
                        )
                        println(io, rpad(line, 80))
                    end
                    println(io, "END")
                end
            end
        end
        return data_root
    end

    system = "toy_system"
    k = 3
    n_conformers = 9
    fixtures = _write_fixture(mktempdir(), system)

    # First pass (no scores) to learn the seeded KMeans labels so the crafted score table is
    # independent of which label maps to which shape group.
    base = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = mktempdir(),
        kmeans_k = k,
        seed = 42,
    )

    # No score table: stage skipped cleanly, pipeline still completes over all conformers.
    @test base.deepaccnet_dir === nothing
    @test base.surviving === nothing
    @test DataFrames.nrow(base.hierarchical) == n_conformers

    label_col = "KMeans_K$k"
    labels = base.clustering[!, label_col]
    names = base.clustering.Name

    # Target the smallest cluster (a strict minority of the ensemble) so that giving its members
    # a low score and everyone else a high score makes only that cluster fail at rel 50 / frac 50.
    target = sort(unique(labels); by = l -> (count(==(l), labels), l))[1]
    target_size = count(==(target), labels)
    survivors = sort([l for l in unique(labels) if l != target])
    survivor_names = Set(names[labels .!= target])

    score_table = DataFrames.DataFrame(
        "Name" => names,
        "cb-lddt" => [lbl == target ? 0.0 : 1.0 for lbl in labels],
    )

    out = mktempdir()
    result = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = out,
        kmeans_k = k,
        seed = 42,
        score_table = score_table,
        surviving_rel = 50,
        surviving_frac = 50,
    )

    # The deepaccnet/ tree holds only the single configured cut — no leftover grid sweep, and the
    # filename no longer encodes the rel/frac values.
    @test result.deepaccnet_dir == joinpath(out, system, "deepaccnet")
    surv_dir = joinpath(result.deepaccnet_dir, "surviving")
    @test isdir(surv_dir)
    @test isfile(joinpath(surv_dir, "surviving.csv"))
    @test readdir(surv_dir) == ["surviving.csv"]

    # The configured cut drops exactly the target cluster.
    @test result.surviving == survivors
    chosen = CSV.read(joinpath(surv_dir, "surviving.csv"), DataFrames.DataFrame)
    @test DataFrames.names(chosen) == ["Name", label_col]
    @test Set(chosen.Name) == survivor_names
    @test all(in(Set(survivors)), chosen[!, label_col])

    # Stage 3 switches to the surviving set: the hierarchical table covers only surviving members.
    @test DataFrames.nrow(result.hierarchical) == n_conformers - target_size
    @test Set(result.hierarchical.Name) == survivor_names
    @test all(in(Set(survivors)), result.hierarchical[!, label_col])

    # A surviving cut outside the valid 0..100 percentile range is a clear error.
    @test_throws ErrorException cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = mktempdir(),
        kmeans_k = k,
        seed = 42,
        score_table = score_table,
        surviving_rel = 150,
        surviving_frac = 50,
    )
end

@testitem "Conformer clustering score filter drops unscored" begin
    import CSV
    import DataFrames
    import Random

    function _write_fixture(data_root, system)
        models = joinpath(data_root, system, "models")
        mkpath(models)
        residues = 10
        base = [(3.8 * (i - 1), 0.0, 0.0) for i = 1:residues]
        for (group, bend) in ["A" => 0.0, "B" => 5.0, "C" => 12.0]
            for rep = 1:3
                rng = Random.MersenneTwister(abs(hash((group, rep))) % 100_000)
                path = joinpath(models, "conf_$(group)$(rep).pdb")
                open(path, "w") do io
                    for (i, b) in enumerate(base)
                        yb =
                            i > residues ÷ 2 ? bend * (i - residues ÷ 2) / (residues ÷ 2) :
                            0.0
                        jitter = 0.05 .* (Random.rand(rng, 3) .- 0.5)
                        x, y, z = b[1] + jitter[1], b[2] + yb + jitter[2], b[3] + jitter[3]
                        line = string(
                            "ATOM  ",
                            lpad(i, 5),
                            "  CA  ALA A",
                            lpad(i, 4),
                            "    ",
                            lpad(string(round(x, digits = 3)), 8),
                            lpad(string(round(y, digits = 3)), 8),
                            lpad(string(round(z, digits = 3)), 8),
                            "  1.00  0.00           C",
                        )
                        println(io, rpad(line, 80))
                    end
                    println(io, "END")
                end
            end
        end
        return data_root
    end

    system = "toy_system"
    k = 3
    fixtures = _write_fixture(mktempdir(), system)

    base = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = mktempdir(),
        kmeans_k = k,
        seed = 42,
    )
    label_col = "KMeans_K$k"
    names = base.clustering.Name

    # All conformers score equally (every cluster survives), but one conformer is left out of the
    # score table. It must be dropped from both the surviving tables and the hierarchical output.
    omitted = names[1]
    kept = names[2:end]
    score_table = DataFrames.DataFrame("Name" => kept, "cb-lddt" => fill(1.0, length(kept)))

    out = mktempdir()
    result = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = out,
        kmeans_k = k,
        seed = 42,
        score_table = score_table,
    )

    @test result.deepaccnet_dir == joinpath(out, system, "deepaccnet")
    surv_csv = joinpath(result.deepaccnet_dir, "surviving", "surviving.csv")
    surviving = CSV.read(surv_csv, DataFrames.DataFrame)
    @test !(omitted in surviving.Name)
    @test Set(surviving.Name) == Set(kept)

    # The hierarchical output covers the scored conformers only; the omitted one is absent.
    @test !(omitted in result.hierarchical.Name)
    @test Set(result.hierarchical.Name) == Set(kept)
end

@testitem "Conformer clustering dimensionality reduction" begin
    # _pca (MultivariateStats): points on a line span one principal axis, so the projection has a
    # single real dimension. Columns are observations.
    line = Float64[-1.5 -0.5 0.5 1.5; -1.5 -0.5 0.5 1.5]   # 2 features, 4 points on y = x
    proj = AlphaConformers._pca(line, 2)
    @test size(proj) == (1, 4)                             # only one real direction (no padding)
    @test maximum(proj[1, :]) - minimum(proj[1, :]) > 0    # all spread on the first axis
    # The projection of evenly spaced points stays evenly spaced (orthogonal transform).
    @test proj[1, 2] - proj[1, 1] ≈ -(proj[1, 3] - proj[1, 4]) atol = 1e-8

    # A single-feature input has a single direction too.
    flat = Float64[1.0 2.0 3.0]                            # 1 feature, 3 points
    proj_flat = AlphaConformers._pca(flat, 2)
    @test size(proj_flat) == (1, 3)

    # _ensure_2d pads a sub-2-D projection to 2×N with a zero row, for plotting.
    padded = AlphaConformers._ensure_2d(proj)
    @test size(padded) == (2, 4)
    @test padded[1, :] == proj[1, :]
    @test all(padded[2, :] .== 0.0)
    already = Float64[1 2; 3 4]
    @test AlphaConformers._ensure_2d(already) === already   # a ≥2-row matrix is returned as-is

    # _classical_mds recovers a known geometry up to rotation/reflection: a 3-4-5 triangle.
    dist = Float64[0 3 4; 3 0 5; 4 5 0]
    coords = AlphaConformers._classical_mds(dist, 2)
    @test size(coords) == (2, 3)
    recovered = [sqrt(sum(abs2, coords[:, i] .- coords[:, j])) for i = 1:3, j = 1:3]
    @test recovered ≈ dist atol = 1e-6

    # An empty distance matrix yields an ndims×0 layout.
    @test AlphaConformers._classical_mds(zeros(0, 0), 2) == zeros(2, 0)

    # Collinear points span a single axis: the second axis collapses to ~0.
    collinear = Float64[0 1 2; 1 0 1; 2 1 0]
    flat_coords = AlphaConformers._classical_mds(collinear, 2)
    @test all(abs.(flat_coords[2, :]) .< 1e-6)

    # _feature_distance_matrix is symmetric with a zero diagonal and exact euclidean entries.
    feats = Float64[0 3 0; 0 0 4]                          # points (0,0), (3,0), (0,4)
    fdist = AlphaConformers._feature_distance_matrix(feats)
    @test fdist == fdist'
    @test all(fdist[i, i] == 0.0 for i = 1:3)
    @test fdist[1, 2] ≈ 3.0
    @test fdist[2, 3] ≈ 5.0

    # _plot_dendrogram with a `nothing` tree (a cluster too small to merge) still writes a file.
    ENV["GKSwstype"] = "100"
    placeholder = joinpath(mktempdir(), "dendrogram.png")
    @test AlphaConformers._plot_dendrogram(nothing, 1.0, "tiny", placeholder) == placeholder
    @test isfile(placeholder)
end

@testitem "Conformer clustering styling and graph helpers" begin
    # Reuse the Plots binding the package already imports, so the test environment needs no
    # extra dependency.
    viridis = AlphaConformers.plt.cgrad(:viridis)

    # _cluster_color_map: the lowest label maps to one gradient end, the highest to the other,
    # at position i/(n-1); a single label maps to the midpoint; a -1 noise label stays gray.
    cmap = AlphaConformers._cluster_color_map([0, 1, 2])
    @test cmap[0] == viridis[0.0]
    @test cmap[1] == viridis[0.5]
    @test cmap[2] == viridis[1.0]
    @test sort(collect(keys(cmap))) == [0, 1, 2]

    single = AlphaConformers._cluster_color_map([5])
    @test single[5] == viridis[0.5]

    noisy = AlphaConformers._cluster_color_map([-1, 0, 1])
    @test noisy[-1] == AlphaConformers._NOISE_COLOR     # light gray
    @test noisy[0] == viridis[0.0]                      # noise is excluded from the gradient
    @test noisy[1] == viridis[1.0]

    # _mst_edges: a path 1-2-3 (weights 1 then 2) is the minimum spanning backbone.
    d = Float64[0 1 5; 1 0 2; 5 2 0]
    @test Set(AlphaConformers._mst_edges(d)) == Set([(1, 2, 1.0), (2, 3, 2.0)])
    @test AlphaConformers._mst_edges(zeros(1, 1)) == Tuple{Int,Int,Float64}[]
    @test AlphaConformers._mst_edges(zeros(0, 0)) == Tuple{Int,Int,Float64}[]

    # _intercluster_rmsd: average-linkage RMSD between sub-clusters (two tight pairs far apart).
    dist = Float64[0 0.5 5 5; 0.5 0 5 5; 5 5 0 0.5; 5 5 0.5 0]
    uniq, mean = AlphaConformers._intercluster_rmsd(dist, [1, 1, 2, 2])
    @test uniq == [1, 2]
    @test mean[1, 1] ≈ 0.25       # within sub-cluster 1: (0 + 0.5 + 0.5 + 0)/4
    @test mean[2, 2] ≈ 0.25
    @test mean[1, 2] ≈ 5.0        # between the two far-apart pairs
    @test mean[1, 2] == mean[2, 1]

    # _cluster_centroid_rmsd: centroid euclidean distance divided by sqrt(n_atoms).
    features = Float64[0 0 3; 0 0 0; 0 0 0]    # 1 atom; two points at origin, one at (3,0,0)
    ids, inter = AlphaConformers._cluster_centroid_rmsd(features, [0, 0, 1], 1)
    @test ids == [0, 1]
    @test inter[1, 2] ≈ 3.0
    @test inter == inter'
    @test all(inter[i, i] == 0.0 for i = 1:2)

    # _plot_minicluster_graph is skipped (returns nothing) when there is only one sub-cluster.
    ENV["GKSwstype"] = "100"
    single_path = joinpath(mktempdir(), "minicluster_graph.png")
    @test AlphaConformers._plot_minicluster_graph(
        zeros(2, 2),
        [1, 1],
        "one sub-cluster",
        single_path,
    ) === nothing
    @test !isfile(single_path)
    # With two sub-clusters it writes the figure.
    multi_path = joinpath(mktempdir(), "minicluster_graph.png")
    @test AlphaConformers._plot_minicluster_graph(dist, [1, 1, 2, 2], "two", multi_path) ==
          multi_path
    @test isfile(multi_path)
end

@testitem "Conformer clustering essential plots" begin
    import DataFrames
    import Random

    # Same deterministic ensemble as the other integration tests (three shape groups).
    function _write_fixture(data_root, system)
        models = joinpath(data_root, system, "models")
        mkpath(models)
        residues = 10
        base = [(3.8 * (i - 1), 0.0, 0.0) for i = 1:residues]
        for (group, bend) in ["A" => 0.0, "B" => 5.0, "C" => 12.0]
            for rep = 1:3
                rng = Random.MersenneTwister(abs(hash((group, rep))) % 100_000)
                path = joinpath(models, "conf_$(group)$(rep).pdb")
                open(path, "w") do io
                    for (i, b) in enumerate(base)
                        yb =
                            i > residues ÷ 2 ? bend * (i - residues ÷ 2) / (residues ÷ 2) :
                            0.0
                        jitter = 0.05 .* (Random.rand(rng, 3) .- 0.5)
                        x, y, z = b[1] + jitter[1], b[2] + yb + jitter[2], b[3] + jitter[3]
                        line = string(
                            "ATOM  ",
                            lpad(i, 5),
                            "  CA  ALA A",
                            lpad(i, 4),
                            "    ",
                            lpad(string(round(x, digits = 3)), 8),
                            lpad(string(round(y, digits = 3)), 8),
                            lpad(string(round(z, digits = 3)), 8),
                            "  1.00  0.00           C",
                        )
                        println(io, rpad(line, 80))
                    end
                    println(io, "END")
                end
            end
        end
        return data_root
    end

    system = "toy_system"
    k = 3
    fixtures = _write_fixture(mktempdir(), system)
    models = joinpath(fixtures, system, "models")
    refs = mktempdir()
    cp(joinpath(models, "conf_A1.pdb"), joinpath(refs, "REFAP_A.pdb"))
    cp(joinpath(models, "conf_C1.pdb"), joinpath(refs, "REFHO_B.pdb"))

    # Figures dropped by the port must never appear under the output tree. (The mini-cluster
    # graph and the cluster-colored query path are re-added in Slice 8, so they are no longer
    # dropped; the 3D PCA, transition-network and benchmark/validation variants stay out.)
    dropped = (
        "elbow",
        "transition",
        "negative_control",
        "spearman",
        "bubble",
        "af3",
        "3d",
        "sweep",
        "benchmark",
        "validation",
    )
    function assert_no_dropped(out)
        for (_root, _dirs, files) in walkdir(out), f in files
            endswith(f, ".png") || continue
            @test !any(d -> occursin(d, lowercase(f)), dropped)
        end
    end

    # Every figure carries the Slice 8 styling (600 dpi, 13×8 inches).
    @test AlphaConformers._FIG_DPI == 600
    @test AlphaConformers._FIG_SIZE == (1300, 800)

    # 0-reference mode: system-level cluster overview + plain embedding scatter, per-cluster
    # agglomerative figures, no reference-only views and no affinity/query views.
    out0 = mktempdir()
    r0 = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = out0,
        kmeans_k = k,
        subcluster_threshold = 0.001,   # split each cluster so the mini-cluster graph is produced
    )
    @test isfile(joinpath(r0.system_dir, "cluster_overview.png"))
    @test isfile(joinpath(r0.system_dir, "reference_scatter.png"))
    @test !isfile(joinpath(r0.system_dir, "reference_affinity.png"))
    @test !isfile(joinpath(r0.system_dir, "query_path.png"))
    # The individual-RMSD and reference-cluster-assignment plots were dropped in this iteration.
    @test !isfile(joinpath(r0.system_dir, "individual_rmsd.png"))
    @test !isfile(joinpath(r0.system_dir, "reference_cluster_assignment.png"))
    for cluster_dir in r0.cluster_dirs
        @test isfile(joinpath(cluster_dir, "subcluster_scatter.png"))
        @test isfile(joinpath(cluster_dir, "dendrogram.png"))
        @test isfile(joinpath(cluster_dir, "minicluster_graph.png"))
        # The RMSD heatmap was dropped in this iteration.
        @test !isfile(joinpath(cluster_dir, "rmsd_heatmap.png"))
    end
    @test all(isfile, r0.figures)                            # the returned paths exist
    assert_no_dropped(out0)

    # 1-reference mode: adaptive reference scatter, still no affinity.
    out1 = mktempdir()
    r1 = cluster_conformers(
        system,
        "REFAP_A";
        refs_dir = refs,
        data_root = fixtures,
        out_dir = out1,
        kmeans_k = k,
    )
    @test isfile(joinpath(r1.system_dir, "cluster_overview.png"))
    @test isfile(joinpath(r1.system_dir, "reference_scatter.png"))
    @test !isfile(joinpath(r1.system_dir, "individual_rmsd.png"))
    @test !isfile(joinpath(r1.system_dir, "reference_affinity.png"))
    # With one reference and no explicit query, the query defaults to that reference, so the
    # query-path view is produced automatically.
    @test isfile(joinpath(r1.system_dir, "query_path.png"))
    @test all(isfile, r1.figures)
    assert_no_dropped(out1)

    # 2-reference mode: adaptive (apo-vs-holo) reference scatter; affinity plot was removed.
    out2 = mktempdir()
    r2 = cluster_conformers(
        system,
        "REFAP_A",
        "REFHO_B";
        refs_dir = refs,
        data_root = fixtures,
        out_dir = out2,
        kmeans_k = k,
    )
    @test isfile(joinpath(r2.system_dir, "cluster_overview.png"))
    @test isfile(joinpath(r2.system_dir, "reference_scatter.png"))
    @test !isfile(joinpath(r2.system_dir, "individual_rmsd.png"))
    @test !isfile(joinpath(r2.system_dir, "reference_affinity.png"))
    # With two references and no explicit query, the query defaults to the first reference, so the
    # query-path view is produced and starts from that reference's cluster.
    @test isfile(joinpath(r2.system_dir, "query_path.png"))
    @test r2.reference_clusters.Reference[1] == "ref1"
    @test all(isfile, r2.figures)
    assert_no_dropped(out2)

    # A supplied query (a conformer name) adds the query-path view (C4).
    out_q = mktempdir()
    rq = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = out_q,
        kmeans_k = k,
        query = r0.clustering.Name[1],
    )
    @test isfile(joinpath(rq.system_dir, "query_path.png"))
    @test any(f -> occursin("query_path", f), rq.figures)

    # A query given as a structure-file path (assigned to the nearest cluster centroid), in
    # 2-reference mode so the reference-assigned clusters are overlaid as stars on the path.
    out_qp = mktempdir()
    rqp = cluster_conformers(
        system,
        "REFAP_A",
        "REFHO_B";
        refs_dir = refs,
        data_root = fixtures,
        out_dir = out_qp,
        kmeans_k = k,
        query = joinpath(models, "conf_B2.pdb"),
    )
    @test isfile(joinpath(rqp.system_dir, "query_path.png"))

    # An unresolvable query is a clear error.
    @test_throws ErrorException cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = mktempdir(),
        kmeans_k = k,
        query = "no_such_conformer",
    )

    # Score-filter mode: the three score-stage figures appear only when the filter runs (C2).
    out_s = mktempdir()
    score_table = DataFrames.DataFrame(
        "Name" => r0.clustering.Name,
        "cb-lddt" => collect(range(0.1, 0.9; length = DataFrames.nrow(r0.clustering))),
    )
    rs = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = out_s,
        kmeans_k = k,
        score_table = score_table,
        surviving_rel = 50,
        surviving_frac = 75,
    )
    @test rs.deepaccnet_dir !== nothing
    @test isfile(joinpath(rs.deepaccnet_dir, "rmsd_scatter_by_score.png"))
    @test isfile(joinpath(rs.deepaccnet_dir, "kmeans_deepaccnet_filtered.png"))
    # The score-distribution histogram was dropped in this iteration.
    @test !isfile(joinpath(rs.deepaccnet_dir, "score_distribution.png"))
    @test all(isfile, rs.figures)
    assert_no_dropped(out_s)
    # The score-stage figures are not produced in a reference-free run with no score table.
    @test !any(f -> occursin("rmsd_scatter_by_score", f), r0.figures)

    # make_plots = false writes no figures at all.
    out_n = mktempdir()
    rn = cluster_conformers(
        system;
        data_root = fixtures,
        out_dir = out_n,
        kmeans_k = k,
        make_plots = false,
    )
    @test isempty(rn.figures)
    @test !isfile(joinpath(rn.system_dir, "cluster_overview.png"))
end

@testitem "Conformer discovery reads only prediction models, not templates" begin
    # The AlphaConformers output tree nests predictions under
    # `cluster_*/.../predictions/sequences/models/` and input templates under
    # `templates_complete/` and `templates_adaptative/`. Discovery must read only the former.
    root = mktempdir()
    system = "sys"
    sys_dir = joinpath(root, system)

    models = joinpath(sys_dir, "cluster_1", "af", "predictions", "sequences", "models")
    mkpath(models)
    touch(joinpath(models, "seed0.pdb"))
    touch(joinpath(models, "seed1.cif"))

    # Templates and search results must be ignored even though they hold .pdb/.cif files.
    for sub in ("templates_complete", "templates_adaptative")
        d = joinpath(sys_dir, "cluster_1", sub)
        mkpath(d)
        touch(joinpath(d, "tmpl.pdb"))
    end
    fullpdb = joinpath(sys_dir, "fullpdb_results")
    mkpath(fullpdb)
    touch(joinpath(fullpdb, "hit.pdb"))

    found = AlphaConformers._discover_conformers(root, system)
    names_found = first.(found)
    @test length(found) == 2
    @test all(occursin(joinpath("models", ""), n) for n in names_found)
    @test !any(occursin("templates", n) for n in names_found)
    @test !any(occursin("fullpdb_results", n) for n in names_found)

    # A system with only templates (no models/) errors clearly instead of clustering them.
    bare = mktempdir()
    tdir = joinpath(bare, system, "templates_complete")
    mkpath(tdir)
    touch(joinpath(tdir, "tmpl.pdb"))
    @test_throws ErrorException AlphaConformers._discover_conformers(bare, system)

    # An explicit `subdir` (with glob wildcards) reads only that method subfolder, skipping the
    # full walk. Both the literal and the wildcard form find exactly the two prediction models.
    literal = AlphaConformers._discover_conformers(
        root,
        system;
        subdir = joinpath("cluster_1", "af", "predictions", "sequences", "models"),
    )
    @test Set(first.(literal)) == Set(names_found)
    wildcard = AlphaConformers._discover_conformers(
        root,
        system;
        subdir = joinpath("cluster_*", "af", "predictions", "sequences", "models"),
    )
    @test Set(first.(wildcard)) == Set(names_found)
    # A subdir that matches nothing errors clearly.
    @test_throws ErrorException AlphaConformers._discover_conformers(
        root,
        system;
        subdir = "no_such_folder",
    )
end

@testitem "cluster_conformers is the documented public entry point" begin
    # The pipeline's only public surface is the exported entry point; the stage
    # functions stay internal (underscore-prefixed, unexported).
    @test :cluster_conformers in names(AlphaConformers)
    @test isdefined(AlphaConformers, :cluster_conformers)

    # The exported binding is the callable pipeline, not a stray alias.
    @test cluster_conformers isa Function

    # It carries a docstring so it is discoverable from the REPL and the docs.
    doc = string(@doc cluster_conformers)
    @test occursin("cluster_conformers", doc)
    @test occursin("# Behavior", doc)

    # The internal stage helpers are not exported.
    @test :_kmeans_cluster ∉ names(AlphaConformers)
    @test :_write_surviving_table ∉ names(AlphaConformers)
end
