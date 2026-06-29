using TestItems

@testitem "USalign command and table helpers" begin
    import DataFrames

    cmd = AlphaConformers._usalign_command(
        `usalign`;
        mol = "RNA",
        TMcut = nothing,
        fast = false,
        outfmt = 1,
        extra_parameters = `-ter 0`,
    )
    cmd_text = string(cmd)
    @test occursin("-mol RNA", cmd_text)
    @test occursin("-outfmt 1", cmd_text)
    @test occursin("-ter 0", cmd_text)
    @test !occursin("-TMcut", cmd_text)
    @test !occursin("-fast", cmd_text)

    default_cmd = string(AlphaConformers._usalign_command(`usalign`))
    @test occursin("-TMcut 0.5", default_cmd)
    @test occursin("-fast", default_cmd)

    @test_throws AssertionError AlphaConformers._usalign_command(`usalign`; mol = "DNA")
    @test_throws AssertionError AlphaConformers._usalign_command(`usalign`; TMcut = 0.1)
    @test_throws AssertionError AlphaConformers._usalign_command(`usalign`; outfmt = 3)

    mktempdir() do dir
        pdb_folder = joinpath(dir, "pdbs")
        mkdir(pdb_folder)
        touch(joinpath(pdb_folder, "target_a.pdb"))
        list_file = joinpath(dir, "list.txt")

        AlphaConformers._save_list(list_file, ["target_a.pdb", "missing.pdb"], pdb_folder)
        @test readlines(list_file) == ["target_a.pdb"]

        table_file = joinpath(dir, "usalign.tsv")
        write(table_file, "#PDBchain1\tPDBchain2\tTM1\nquery\t target\t1.0\n")
        df = AlphaConformers._read_usalign_output_table(table_file)
        @test DataFrames.nrow(df) == 1
        @test "PDBchain1" in names(df)
        @test df[1, :PDBchain1] == "query"

        query = joinpath(dir, "query.pdb")
        touch(query)
        target = joinpath(dir, "target.pdb")
        touch(target)
        fake_usalign = joinpath(dir, "fake_usalign")
        write(
            fake_usalign,
            "#!/bin/sh\n" *
            "printf '#PDBchain1\tPDBchain2\tTM1\\nquery\\ttarget\\t1.0\\n'\n",
        )
        chmod(fake_usalign, 0o755)
        fake_cmd = Cmd([fake_usalign])

        pairwise = usalign(query, target; usalign = fake_cmd)
        @test DataFrames.nrow(pairwise) == 1
        @test pairwise[1, :PDBchain1] == "query"

        batch = usalign(query, pdb_folder, ["target_a.pdb"]; usalign = fake_cmd)
        @test DataFrames.nrow(batch) == 1

        @test_throws AssertionError usalign(joinpath(dir, "missing.pdb"), query)
        @test_throws AssertionError usalign(query, joinpath(dir, "missing.pdb"))
        @test_throws AssertionError usalign(
            query,
            joinpath(dir, "missing_folder"),
            ["x.pdb"],
        )
    end
end
