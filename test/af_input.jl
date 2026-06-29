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

@testitem "AlphaFold input readers cover mmCIF, lowercase chains, and error paths" begin
    import MIToS
    import Printf

    function write_test_pdb(path; chain = "A")
        residues = ["ALA", "GLY", "SER"]
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
        lower_pdb = joinpath(dir, "lower.pdb")
        write_test_pdb(lower_pdb; chain = "a")

        lower_chain = AlphaConformers._read_pdb_chain(lower_pdb, "A")
        @test lower_chain !== nothing
        @test only(unique(r.id.chain for r in lower_chain)) == "a"
        @test AlphaConformers._read_pdb_chain(lower_pdb, "Z") === nothing
        @test AlphaConformers._read_pdb_chain(joinpath(dir, "missing.pdb")) === nothing
        @test AlphaConformers._read_pdb_chain(joinpath(dir, "missing.cif")) === nothing
        @test AlphaConformers._read_pdb_chain(joinpath(dir, "missing.pdb"), "A") === nothing
        @test AlphaConformers._read_pdb_chain(joinpath(dir, "missing.cif"), "A") === nothing

        cif = joinpath(dir, "lower.cif")
        all_residues = AlphaConformers._read_pdb_chain(lower_pdb)
        MIToS.PDB.write_file(cif, all_residues, MIToS.PDB.MMCIFFile)
        @test AlphaConformers._read_pdb_chain(cif) !== nothing
        @test AlphaConformers._read_pdb_chain(cif, "A") !== nothing

        data = get_residues_and_sequence(lower_pdb; chain = "A")
        @test data.sequence == "AGS"
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

@testitem "Sequence and template input helpers write expected lightweight files" begin
    import MIToS
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
        ref = joinpath(dir, "ref.pdb")
        template = joinpath(dir, "template.pdb")
        write_test_pdb(ref)
        write_test_pdb(template)

        paths = create_pdb_lists(ref, "A", "1", [template, ref], ["A", "A"], ["1", "1"])
        @test paths.pdb_files[1] == abspath(ref)
        @test paths.pdb_files[2] == abspath(template)

        paths_missing_ref = create_pdb_lists(ref, "A", "1", [template], ["A"], ["1"])
        @test paths_missing_ref.pdb_files[1] == abspath(ref)
        @test paths_missing_ref.chains[1] == "A"

        save_sequences(dir, [ref, template]; chains = ["A", "A"], models = ["1", "1"])
        fasta = read(joinpath(dir, "sequences.fasta"), String)
        @test occursin(">ref.pdb", fasta)
        @test occursin(">template.pdb", fasta)
        @test occursin("AGST", fasta)

        msa_file = joinpath(dir, "small.fasta")
        write(msa_file, ">query\nAGST\n>template\nAG-T\n")
        msa = read(msa_file, MIToS.MSA.FASTA)
        cleaned = clean_msa(msa)
        @test size(cleaned, 1) == 2
        @test size(cleaned, 2) == 4

        aligned = align_sequences(msa_file)
        @test size(aligned, 1) == 2
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
