using TestItems

@testitem "PDB folder helpers use local databases and dispatch target containers" begin
    import DataFrames
    import Printf

    function write_test_pdb(path; chain = "A")
        residues = ["ALA", "GLY", "SER", "THR"]
        serial = 1
        open(path, "w") do io
            for (resnum, resname) in enumerate(residues)
                x = 3.8 * (resnum - 1)
                atoms = (
                    ("N", x - 1.2, 0.0, 0.0, "N"),
                    ("CA", x, 0.0, 0.0, "C"),
                    ("C", x + 1.2, 1.0, 0.0, "C"),
                    ("O", x + 1.8, 2.0, 0.0, "O"),
                )
                for (atom, ax, ay, az, element) in atoms
                    Printf.@printf(
                        io,
                        "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n",
                        serial,
                        atom,
                        resname,
                        chain,
                        resnum,
                        ax,
                        ay,
                        az,
                        element,
                    )
                    serial += 1
                end
            end
            println(io, "TER")
            println(io, "END")
        end
    end

    mktempdir() do dir
        pdb_db = joinpath(dir, "pdb_db")
        mkdir(pdb_db)
        write_test_pdb(joinpath(pdb_db, "1ABC.pdb"))
        write_test_pdb(joinpath(pdb_db, "2DEF.pdb"))
        write_test_pdb(joinpath(pdb_db, "AF-P11111-F1-model_v4.pdb"))

        chain = AlphaConformers._get_pdb_chain(pdb_db, "1ABC", "A")
        @test chain !== nothing
        @test length(chain) > 0
        @test AlphaConformers._get_pdb_chain(pdb_db, "9ZZZ", "A") === nothing
        @test AlphaConformers._get_pdb_chain(
            "",
            "1ABC",
            "A";
            download_pdb_chain = (pdb, chain) -> (pdb, chain),
        ) == ("1ABC", "A")

        downloaded_path = joinpath(dir, "downloaded.pdb.gz")
        fake_download(pdb; filename, format) = write(filename, "downloaded $pdb")
        function fake_read(path, chain)
            @test path == downloaded_path
            @test chain == "A"
            @test isfile(path)
            return (path, chain)
        end
        @test AlphaConformers._download_pdb_chain(
            "1ABC",
            "A";
            downloadpdb = fake_download,
            read_pdb_chain = fake_read,
            tmp_path = downloaded_path,
        ) == (downloaded_path, "A")
        @test !isfile(downloaded_path)

        partial_path = joinpath(dir, "partial.pdb.gz")
        function failing_download(pdb; filename, format)
            write(filename, "partial $pdb")
            error("download failed")
        end
        @test AlphaConformers._download_pdb_chain(
            "1ABC",
            "A";
            downloadpdb = failing_download,
            read_pdb_chain = fake_read,
            tmp_path = partial_path,
        ) === nothing
        @test !isfile(partial_path)

        array_dir = create_pdb_folder(["1ABC.pdb_A"], joinpath(dir, "array"); pdb_db)
        @test isfile(joinpath(array_dir, "1ABC.pdb_A"))

        df = DataFrames.DataFrame(target = ["2DEF.pdb_A"])
        df_dir = create_pdb_folder(df, joinpath(dir, "dataframe"); pdb_db)
        @test isfile(joinpath(df_dir, "2DEF.pdb_A"))

        af_dir = create_pdb_folder(
            Set(["AF-P11111-F1-model_v4.pdb"]),
            joinpath(dir, "af");
            pdb_db,
        )
        @test isfile(joinpath(af_dir, "AF-P11111-F1-model_v4.pdb"))

        missing_dir = create_pdb_folder(["9ZZZ.pdb_A"], joinpath(dir, "missing"); pdb_db)
        @test isdir(missing_dir)
        @test !isfile(joinpath(missing_dir, "9ZZZ.pdb_A"))
    end
end

@testitem "Empty folder helper replaces existing contents" begin
    mktempdir() do dir
        folder = joinpath(dir, "scratch")
        mkdir(folder)
        write(joinpath(folder, "old.txt"), "old")

        AlphaConformers._create_empty_folder(folder)

        @test isdir(folder)
        @test isempty(readdir(folder))

        nested = joinpath(dir, "missing", "nested", "scratch")
        AlphaConformers._create_empty_folder(nested)

        @test isdir(nested)
        @test isempty(readdir(nested))
    end
end
