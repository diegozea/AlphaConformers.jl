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

    @test AlphaConformers._container_runtime("apptainer") == "apptainer"
    @test AlphaConformers._container_runtime("singularity") == "singularity"
    @test_throws ArgumentError AlphaConformers._container_runtime("")
    @test_throws ArgumentError AlphaConformers._container_runtime(" apptainer ")
    @test_throws ArgumentError AlphaConformers._container_runtime("/opt/site/bin/container")

    cmd = AlphaConformers._colabfold_batch_command(
        "/tmp/cluster_1",
        "/tmp/cluster_1/af",
        "/tmp/colabfold.sif",
        "/tmp/colabfold-cache";
        seed = 12_345,
        tmp_dir = "/tmp/cluster_1/af/tmp",
    )
    @test cmd.exec[1] == "apptainer"
    @test cmd.exec[2] == "exec"
    @test "/tmp/colabfold-cache:/cache/colabfold" in cmd.exec
    @test "/tmp/cluster_1/af/tmp:/mnt/tmp" in cmd.exec
    @test any(arg -> occursin("TMPDIR=/mnt/tmp", arg), cmd.exec)
    @test "colabfold_batch" in cmd.exec

    singularity_cmd = AlphaConformers._colabfold_batch_command(
        "/tmp/cluster_1",
        "/tmp/cluster_1/af",
        "/tmp/colabfold.sif",
        "/tmp/colabfold-cache";
        seed = 12_345,
        tmp_dir = "/tmp/cluster_1/af/tmp",
        container_runtime = "singularity",
    )
    @test singularity_cmd.exec[1] == "singularity"
    @test singularity_cmd.exec[2] == "exec"

    @test_throws ArgumentError AlphaConformers._colabfold_batch_command(
        "/tmp/cluster_1",
        "/tmp/cluster_1/af",
        "/tmp/colabfold.sif",
        "/tmp/colabfold-cache";
        seed = 12_345,
        tmp_dir = "/tmp/cluster_1/af/tmp",
        container_runtime = "",
    )

    @test_throws ArgumentError AlphaConformers._colabfold_batch_command(
        "/tmp/cluster_1",
        "/tmp/cluster_1/af",
        "/tmp/colabfold.sif",
        "/tmp/colabfold-cache";
        seed = 12_345,
        tmp_dir = "/tmp/cluster_1/af/tmp",
        container_runtime = "/opt/site/bin/container",
    )

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
