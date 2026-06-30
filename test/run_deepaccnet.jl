using TestItems

# DeepAccNet Apptainer runner
# ---------------------------
# Spec/contract for `run_deepaccnet` and its internal helpers (Slice 9). CI has no GPU or
# apptainer, so — exactly as `_colabfold_batch_command` is tested without executing — every
# part is covered except the handful of lines gated behind an actual container run
# (`run_cmd(cmd; check=true)` and the public `run_deepaccnet` glue that wraps it in
# `mktempdir`): helper-built symlink pairs, the apptainer command string, the dry-run path,
# the input-validation errors, method auto-detection, and symlink creation.

@testitem "run_deepaccnet helpers and orchestration" begin
    # Build one cluster_<N>/<inner>/predictions/sequences/models/ tree with the given PDB
    # file names. Returns the models dir so the test can assert the exact symlink targets.
    function _make_models(base_dir, cluster, inner, pdbs)
        models = joinpath(base_dir, cluster, inner, "predictions", "sequences", "models")
        mkpath(models)
        for f in pdbs
            write(joinpath(models, f), "")
        end
        return models
    end

    # Link-pair construction: naming, alphabetical ordering, filtering, method auto-detection.
    # ------------------------------------------------------------------------------------------
    mktempdir() do data_root
        system = "sys1"
        base = joinpath(data_root, system)
        inner = "af_pdb"
        m1 = _make_models(base, "cluster_1", inner, ["m1.pdb", "m2.pdb"])
        m2 = _make_models(base, "cluster_2", inner, ["m1.pdb"])
        m10 = _make_models(base, "cluster_10", inner, ["m1.pdb"])
        # A non-pdb decoy in a real cluster, and a cluster with no models dir, are ignored.
        write(joinpath(m1, "notes.txt"), "")
        mkpath(joinpath(base, "cluster_3", inner, "predictions", "sequences"))
        # A non-cluster directory is ignored.
        mkpath(joinpath(base, "logs"))

        pairs = AlphaConformers._deepaccnet_links(data_root, system)
        names = first.(pairs)
        targets = last.(pairs)

        # Clusters are matched and visited by name, in plain alphabetical order
        # ("cluster_1" < "cluster_10" < "cluster_2").
        @test names == [
            "cluster_1_m1.pdb",
            "cluster_1_m2.pdb",
            "cluster_10_m1.pdb",
            "cluster_2_m1.pdb",
        ]
        @test targets[1] == joinpath(m1, "m1.pdb")
        @test targets[3] == joinpath(m10, "m1.pdb")
        @test targets[4] == joinpath(m2, "m1.pdb")
        @test all(endswith(".pdb"), names)
        @test !any(n -> occursin("notes.txt", n), names)
    end

    # Method auto-detection refuses an ambiguous system with .pdb files under two different
    # inner method folders.
    # ---------------------------------------------------------------------------------------
    mktempdir() do data_root
        system = "sys_ambiguous"
        base = joinpath(data_root, system)
        _make_models(base, "cluster_1", "af_pdb", ["m1.pdb"])
        _make_models(base, "cluster_1", "af_cif_as_pdb", ["m1.pdb"])

        @test_throws ArgumentError AlphaConformers._deepaccnet_links(data_root, system)
    end

    # `_discover_conformers`'s own "nothing found at all" error is caught and rethrown as an
    # ArgumentError when a cluster_* directory exists but holds no prediction structures
    # whatsoever (not even a `models/` folder).
    # ------------------------------------------------------------------------------------------
    mktempdir() do data_root
        system = "sys_no_predictions"
        mkpath(joinpath(data_root, system, "cluster_1"))

        @test_throws ArgumentError AlphaConformers._deepaccnet_links(data_root, system)
    end

    # Symlink creation: stale paths at the link location are replaced.
    # -----------------------------------------------------------------
    mktempdir() do links_dir
        target1 = joinpath(links_dir, "target1.pdb")
        target2 = joinpath(links_dir, "target2.pdb")
        write(target1, "one")
        write(target2, "two")
        link_path = joinpath(links_dir, "cluster_1_m1.pdb")
        write(link_path, "stale")

        AlphaConformers._symlink_deepaccnet_pairs(
            links_dir,
            [("cluster_1_m1.pdb", target1), ("cluster_1_m2.pdb", target2)],
        )

        @test islink(link_path)
        @test readlink(link_path) == target1
        @test islink(joinpath(links_dir, "cluster_1_m2.pdb"))
    end

    # Apptainer command construction (no execution).
    # ----------------------------------------------
    cmd = AlphaConformers._deepaccnet_command(
        "/tmp/deepaccnet.sif",
        "/data/system",
        "/tmp/run/links",
        "/tmp/run/out/deepaccnet_results.csv";
        n_threads = 4,
    )
    @test "--nv" in cmd.exec
    @test "apptainer" in cmd.exec
    @test "/tmp/deepaccnet.sif" in cmd.exec
    # Data dirs bound: the base dir holding link targets, the symlink dir, and the output dir.
    @test "/data/system" in cmd.exec                   # base dir bound (symlink targets resolve)
    @test "/tmp/run/links" in cmd.exec                 # bound + passed positionally as input
    @test "/tmp/run/out" in cmd.exec                   # output dir bound (dirname of the csv)
    @test "/tmp/run/out/deepaccnet_results.csv" in cmd.exec  # passed positionally as output
    @test count(==("--bind"), cmd.exec) >= 3
    # Featurization CPU count forwarded to the image (`--process N`, after the positionals).
    @test "--process" in cmd.exec
    @test "4" in cmd.exec

    # container_runtime is validated and forwarded as the invoked executable.
    cmd_singularity = AlphaConformers._deepaccnet_command(
        "/tmp/deepaccnet.sif",
        "/data/system",
        "/tmp/run/links",
        "/tmp/run/out/deepaccnet_results.csv";
        container_runtime = "singularity",
    )
    @test "singularity" in cmd_singularity.exec
    @test_throws ArgumentError AlphaConformers._deepaccnet_command(
        "/tmp/deepaccnet.sif",
        "/data/system",
        "/tmp/run/links",
        "/tmp/run/out/deepaccnet_results.csv";
        container_runtime = "docker",
    )

    # Dry run, sample limiting: prints the plan, runs and creates nothing.
    # ---------------------------------------------------------------------
    mktempdir() do links_dir
        pairs = [
            ("cluster_1_m1.pdb", "/data/system/cluster_1/af/m1.pdb"),
            ("cluster_1_m2.pdb", "/data/system/cluster_1/af/m2.pdb"),
        ]
        AlphaConformers._run_deepaccnet(
            links_dir,
            "/data/system",
            "/tmp/deepaccnet.sif",
            joinpath(links_dir, "deepaccnet_results.csv"),
            pairs;
            test_samples = 1,
            dry_run = true,
        )
        @test isempty(readdir(links_dir))
    end

    # Input validation: clear ArgumentErrors before any container work.
    # ----------------------------------------------------------------
    mktempdir() do data_root
        @test_throws ArgumentError AlphaConformers.run_deepaccnet(
            "no_such_system",
            "/tmp/deepaccnet.sif";
            data_root = data_root,
        )
    end
    mktempdir() do data_root
        system = "sys_empty"
        # No cluster_* directories at all.
        mkpath(joinpath(data_root, system))
        @test_throws ArgumentError AlphaConformers.run_deepaccnet(
            system,
            "/tmp/deepaccnet.sif";
            data_root = data_root,
        )
    end
    mktempdir() do data_root
        system = "sys_no_pdb"
        # cluster_* exists, but no *.pdb anywhere under it.
        base = joinpath(data_root, system)
        mkpath(joinpath(base, "cluster_1", "af_cif", "predictions", "sequences", "models"))
        write(
            joinpath(
                base,
                "cluster_1",
                "af_cif",
                "predictions",
                "sequences",
                "models",
                "m1.cif",
            ),
            "",
        )
        @test_throws ArgumentError AlphaConformers.run_deepaccnet(
            system,
            "/tmp/deepaccnet.sif";
            data_root = data_root,
        )
    end
end
