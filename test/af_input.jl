using TestItems

@testitem "AlphaFold input PDB reading helpers" begin
    import MIToS
    import Printf

    function write_test_pdb(path)
        residues = ["ALA", "GLY", "SER", "THR", "LEU"]
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
                        "ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n",
                        serial,
                        atom,
                        resname,
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
        pdb_file = joinpath(dir, "1ABC.pdb")
        write_test_pdb(pdb_file)

        chain_data = get_residues_and_sequence(pdb_file; chain = "A")
        @test !isempty(chain_data.residues)
        @test !isempty(chain_data.sequence)

        all_data = get_residues_and_sequence(pdb_file; chain = MIToS.PDB.All)
        @test !isempty(all_data.residues)
        @test all_data.sequence == chain_data.sequence
    end
end

@testitem "create_pdb_folder reads local PDB targets" begin
    import Printf

    function write_test_pdb(path)
        residues = ["ALA", "GLY", "SER", "THR", "LEU"]
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
                        "ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n",
                        serial,
                        atom,
                        resname,
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

        output_dir = create_pdb_folder(["1ABC.pdb_A"], joinpath(dir, "output"); pdb_db)
        output_file = joinpath(output_dir, "1ABC.pdb_A")

        @test isfile(output_file)
        @test filesize(output_file) > 0
    end
end

@testitem "AlphaFold input query filter" begin
    import DataFrames

    mapping = DataFrames.DataFrame(
        PDB = ["1abc", "1abc", "2def"],
        CHAIN = ["A", "B", "A"],
        SP_PRIMARY = ["P11111", "P11111", "P22222"],
    )
    foldseek_results =
        DataFrames.DataFrame(target = ["1ABC.cif_A", "1ABC.cif_B", "2DEF.cif_A"])

    filtered = AlphaConformers.delete_query_from_target(
        deepcopy(foldseek_results),
        mapping,
        "1ABC",
        "A",
    )

    @test filtered.target == ["2DEF.cif_A"]
    @test DataFrames.nrow(foldseek_results) == 3
end
