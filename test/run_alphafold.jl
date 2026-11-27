using TestItems

@testitem "AlphaFold helper functions" begin
    import JSON3

    @test AlphaConformers.aligned_indices("A-CDE", "ABC-E") == ([0, 1, 3], [0, 2, 3])

    parsed =
        AlphaConformers.parse_plddt_info_in_line("rank_001_model_1 pLDDT=91.23 pTM=0.76")
    @test parsed.rank == "rank_001_model_1"
    @test parsed.pLDDT == 91.23
    @test parsed.pTM == 0.76
    @test AlphaConformers.parse_plddt_info_in_line("no rank here") === nothing

    @test AlphaConformers.run_cmd(`true`) === nothing
    @test AlphaConformers.run_cmd(`false`) === nothing

    @test_throws ErrorException AlphaConformers.run_cmd(`false`; check = true)

    cmd = AlphaConformers._colabfold_batch_command(
        "/tmp/cluster_1",
        "/tmp/cluster_1/af",
        "/tmp/colabfold.sif",
        "/tmp/colabfold-cache";
        seed = 12_345,
        tmp_dir = "/tmp/cluster_1/af/tmp",
    )
    @test "/tmp/colabfold-cache:/cache/colabfold" in cmd.exec
    @test "/tmp/cluster_1/af/tmp:/mnt/tmp" in cmd.exec
    @test any(arg -> occursin("TMPDIR=/mnt/tmp", arg), cmd.exec)
    @test "colabfold_batch" in cmd.exec

    mktempdir() do dir
        @test_throws ArgumentError AlphaConformers.run_alphafold(
            dir,
            "/tmp/colabfold.sif",
            "/tmp/colabfold-cache",
        )
    end

    mktempdir() do dir
        out_json = joinpath(dir, "af3.json")
        result = AlphaConformers.write_af3_json(
            out_json;
            run_name = "query",
            sequence = "ACDE",
            chain_id = "B",
            templates_info = Any[],
            seed = 10,
        )
        @test result == out_json
        config = JSON3.read(read(out_json, String))
        @test config.name == "query"
        @test collect(config.modelSeeds) == [10, 11, 12, 13, 14]
        @test config.sequences[1].protein.id == "B"
        @test config.sequences[1].protein.sequence == "ACDE"
        @test config.dialect == "alphafold3"
    end
end

@testitem "prediction runners preserve existing-output skip guards" begin
    mktempdir() do dir
        # run_alphafold_input_structure skips when af_input/predictions/sequences/models exists
        af_input_root = joinpath(dir, "af_input_run")
        mkpath(joinpath(af_input_root, "af_input", "predictions", "sequences", "models"))
        msa_src = joinpath(dir, "src.a3m")
        write(msa_src, ">query\nACD\n")
        res = redirect_stdout(devnull) do
            AlphaConformers.run_alphafold_input_structure(
                "P12345",
                af_input_root,
                "fake.sif",
                "fake_cache";
                msa_path_save = msa_src,
            )
        end
        @test res === nothing

        # run_alphafold3 skips when cluster_1/af3/predictions/sequences/models exists
        af3_root = joinpath(dir, "1ABC_af3")
        mkpath(joinpath(af3_root, "cluster_1", "af3", "predictions", "sequences", "models"))
        @test redirect_stdout(devnull) do
            AlphaConformers.run_alphafold3(
                af3_root,
                "fake_sif_dir",
                "fake_models",
                "fake_db",
            )
        end === nothing

        # run_boltz2 skips when cluster_1/bz/predictions/sequences/models is non-empty
        bz_root = joinpath(dir, "1ABC_bz")
        c1 = joinpath(bz_root, "cluster_1")
        mkpath(joinpath(c1, "bz", "predictions", "sequences", "models"))
        write(joinpath(c1, "sequences.a3m"), ">1ABC_A\nACD\n")
        write(
            joinpath(c1, "bz", "predictions", "sequences", "models", "model_1.cif"),
            "data_model\n",
        )
        @test redirect_stdout(devnull) do
            AlphaConformers.run_boltz2(bz_root, "fake.sif", "fake_cache")
        end === nothing
        @test isfile(joinpath(c1, "boltz_config.yaml"))
    end
end
