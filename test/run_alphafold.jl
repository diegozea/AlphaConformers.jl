using TestItems

@testitem "AlphaFold helper functions" begin
    import JSON3

    @test AlphaConformers.aligned_indices("A-CDE", "ABC-E") == ([0, 1, 3], [0, 2, 3])

    parsed =
        AlphaConformers.parse_plddt_info_in_line("rank_001_model_1 pLDDT=91.23 pTM=0.76")
    @test parsed.rank == "rank_001_model_1"
    @test parsed.pLDDT == 91.23
    @test parsed.pTM == 0.76

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
