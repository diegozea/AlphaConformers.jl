using TestItems

@testitem "SIFTS mapping helpers" begin
    import DataFrames
    import MIToS

    mapping = DataFrames.DataFrame(
        PDB = ["1abc", "1abc", "2def", "9zzz"],
        CHAIN = ["A", "B", "A", "A"],
        SP_PRIMARY = ["P11111", "P11111", "P22222", "P99999"],
    )

    @test get_uniprot_acc(mapping, "1ABC", "A") == ["P11111"]
    @test get_uniprot_acc(mapping, "1abc") == ["P11111"]
    @test get_uniprot_acc(mapping, "AF-P11111-F1-model_v4.pdb") === nothing

    pdbs = get_pdb_codes(mapping, "P11111")
    @test names(pdbs) == ["PDB", "CHAIN"]
    @test Set(zip(pdbs.PDB, pdbs.CHAIN)) == Set([("1ABC", "A"), ("1ABC", "B")])

    @test AlphaConformers._get_pdb_and_chain("4f4j.pdb_A") == ("4F4J", "A")
    @test AlphaConformers._get_pdb_and_chain("4f4j.pdb") == ("4F4J", MIToS.PDB.All)
    @test AlphaConformers._is_chain("4f4j.pdb_A", "4F4J", "A")

    mktempdir() do dir
        output_path = joinpath(dir, "1AKZ.cif")
        download_call = Ref{Tuple{String,String}}()
        function fake_download(url, path)
            download_call[] = (url, path)
            write(path, "downloaded")
        end

        AlphaConformers._download_rcsb_mmcif(
            "1AKZ.cif",
            output_path;
            download_file = fake_download,
        )

        @test download_call[] == ("https://files.rcsb.org/download/1AKZ.cif", output_path)
        @test isfile(output_path)
    end
end

@testitem "Delete query conformations from Foldseek targets" begin
    import DataFrames

    mapping = DataFrames.DataFrame(
        PDB = ["1abc", "1abc", "2def"],
        CHAIN = ["A", "B", "A"],
        SP_PRIMARY = ["P11111", "P11111", "P22222"],
    )
    search_results = DataFrames.DataFrame(
        target = [
            "1ABC.cif_A",
            "1ABC.cif_B",
            "2DEF.cif_A",
            "AF-P11111-F1-model_v4.pdb",
            "9ZZZ.cif_A",
        ],
    )

    filtered = AlphaConformers.delete_query_from_target(
        deepcopy(search_results),
        mapping,
        "1ABC",
        "A",
    )

    @test filtered.target == ["2DEF.cif_A", "9ZZZ.cif_A"]
    @test DataFrames.nrow(search_results) == 5
end
