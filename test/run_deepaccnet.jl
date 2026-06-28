using TestItems

# DeepAccNet Apptainer runner
# ---------------------------
# Spec/contract for `run_deepaccnet` and its internal helpers (Slice 9). CI has no GPU or
# apptainer, so — exactly as `_colabfold_batch_command` is tested without executing — every
# part is covered except the single GPU-only `run_cmd(cmd; check=true)` line: helper-built
# symlink pairs, the apptainer command string, the dry-run path, the input-validation errors,
# and the intermediates cleanup.

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

    # Link-pair construction: naming, alphabetical ordering, filtering.
    # ------------------------------------------------------------------
    mktempdir() do base
        inner = "af_pdb"
        m1 = _make_models(base, "cluster_1", inner, ["m1.pdb", "m2.pdb"])
        m2 = _make_models(base, "cluster_2", inner, ["m1.pdb"])
        m10 = _make_models(base, "cluster_10", inner, ["m1.pdb"])
        # A non-pdb decoy in a real cluster, and a cluster with no models dir, are ignored.
        write(joinpath(m1, "notes.txt"), "")
        mkpath(joinpath(base, "cluster_3", inner, "predictions", "sequences"))
        # A non-cluster directory is ignored.
        mkpath(joinpath(base, "logs"))

        pairs = AlphaConformers._deepaccnet_links(base, inner)
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

    # Apptainer command construction (no execution).
    # ----------------------------------------------
    cmd = AlphaConformers._deepaccnet_command(
        "/tmp/deepaccnet.sif",
        "/data/system",
        "/tmp/run/links",
        "/tmp/run/out/deepaccnet_results.csv";
        process = 4,
    )
    @test "--nv" in cmd.exec
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

    # Dry run: returns the output csv path, runs nothing, writes nothing.
    # ------------------------------------------------------------------
    mktempdir() do base
        inner = "af"
        _make_models(base, "cluster_1", inner, ["m1.pdb"])
        out_dir = joinpath(base, "deepaccnet_output")

        csv = AlphaConformers.run_deepaccnet(
            base,
            "/tmp/deepaccnet.sif";
            inner = inner,
            output_dir = out_dir,
            dry_run = true,
        )
        @test csv == joinpath(out_dir, "deepaccnet_results.csv")
        # Dry run must not create the links dir, the output dir, or the result csv.
        @test !isdir(joinpath(out_dir, "links"))
        @test !ispath(csv)
    end

    # GPU gate: `_gpu_available` reports a Bool, and a real (non-dry) run refuses up front
    # when no GPU is present (DeepAccNet inference needs `apptainer run --nv`). On a no-GPU
    # host (typical CI) this exercises the refusal; on a GPU host the refusal is skipped (the
    # run would then proceed to the container, which we do not execute here). `dry_run`
    # bypasses the gate so the command/link plan can be inspected anywhere.
    # ----------------------------------------------------------------------------------------
    @test AlphaConformers._gpu_available() isa Bool
    mktempdir() do base
        inner = "af_pdb"
        _make_models(base, "cluster_1", inner, ["m1.pdb"])
        out_dir = joinpath(base, "deepaccnet_output")
        if !AlphaConformers._gpu_available()
            @test_throws ArgumentError AlphaConformers.run_deepaccnet(
                base,
                "/tmp/deepaccnet.sif";
                inner = inner,
                output_dir = out_dir,
            )
            # The refused run must not have created links or output.
            @test !isdir(joinpath(out_dir, "links"))
        end
        # `require_gpu = false` overrides the gate; with `dry_run` it returns without running.
        csv = AlphaConformers.run_deepaccnet(
            base,
            "/tmp/deepaccnet.sif";
            inner = inner,
            output_dir = out_dir,
            require_gpu = false,
            dry_run = true,
        )
        @test csv == joinpath(out_dir, "deepaccnet_results.csv")
    end

    # Input validation: clear ArgumentErrors before any container work.
    # ----------------------------------------------------------------
    @test_throws ArgumentError AlphaConformers.run_deepaccnet(
        "/no/such/base/dir",
        "/tmp/deepaccnet.sif",
    )
    mktempdir() do base
        # No cluster_* directories at all.
        @test_throws ArgumentError AlphaConformers.run_deepaccnet(
            base,
            "/tmp/deepaccnet.sif";
            dry_run = true,
        )
    end
    mktempdir() do base
        # cluster_* exists, but no *.pdb for the requested inner.
        _make_models(base, "cluster_1", "af_pdb", ["m1.pdb"])
        @test_throws ArgumentError AlphaConformers.run_deepaccnet(
            base,
            "/tmp/deepaccnet.sif";
            inner = "af_cif",
            dry_run = true,
        )
    end

    # Intermediates cleanup: heavy files go, results/symlinks stay.
    # ------------------------------------------------------------
    mktempdir() do links
        write(joinpath(links, "bert_x.npy"), "")
        write(joinpath(links, "a.features.npz"), "")
        write(joinpath(links, "b.fa"), "")
        write(joinpath(links, "keep.pdb"), "")

        removed = AlphaConformers._clean_deepaccnet_intermediates(links)

        @test removed == 3
        @test !ispath(joinpath(links, "bert_x.npy"))
        @test !ispath(joinpath(links, "a.features.npz"))
        @test !ispath(joinpath(links, "b.fa"))
        @test ispath(joinpath(links, "keep.pdb"))
    end
end
