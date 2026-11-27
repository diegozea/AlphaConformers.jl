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
        download_call = Ref{Any}()
        function fake_downloadpdb(pdb_code; filename, format)
            download_call[] = (pdb_code, filename, format)
            downloaded_path = endswith(filename, ".gz") ? filename : filename * ".gz"
            write(downloaded_path, "downloaded")
            downloaded_path
        end

        downloaded_path = AlphaConformers._download_rcsb_mmcif(
            "1AKZ.cif",
            output_path;
            downloadpdb = fake_downloadpdb,
        )

        @test download_call[] == ("1AKZ", output_path, MIToS.PDB.MMCIFFile)
        @test downloaded_path == output_path * ".gz"
        @test isfile(downloaded_path)

        gz_output_path = joinpath(dir, "1AKZ.cif.gz")
        gz_downloaded_path = AlphaConformers._download_rcsb_mmcif(
            "1akz.cif.gz",
            gz_output_path;
            downloadpdb = fake_downloadpdb,
        )

        @test download_call[] == ("1AKZ", gz_output_path, MIToS.PDB.MMCIFFile)
        @test gz_downloaded_path == gz_output_path
        @test isfile(gz_downloaded_path)
    end

    mktempdir() do dir
        raw_downloads_dir = joinpath(dir, "raw_downloads")
        tmp_targets_dir = joinpath(dir, "tmp_targets_dir")
        mkdir(tmp_targets_dir)
        download_call = Ref{Any}()

        function fake_download_rcsb_mmcif(prot_name, output_path)
            download_call[] = (prot_name, output_path)
            downloaded_path = output_path * ".gz"
            write(downloaded_path, "downloaded")
            downloaded_path
        end

        downloaded_path = AlphaConformers._download_raw_rcsb_mmcif(
            "1AKZ.cif",
            raw_downloads_dir;
            download_rcsb_mmcif = fake_download_rcsb_mmcif,
        )

        @test download_call[] == ("1AKZ.cif", joinpath(raw_downloads_dir, "1AKZ.cif"))
        @test downloaded_path == joinpath(raw_downloads_dir, "1AKZ.cif.gz")
        @test isfile(downloaded_path)
        @test isempty(readdir(tmp_targets_dir))
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
