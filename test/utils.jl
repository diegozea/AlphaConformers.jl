using TestItems

@testitem "Utility parsers and A3M reader handle edge cases" begin
    mktempdir() do dir
        a3m = joinpath(dir, "test.a3m")
        write(a3m, ">query\nACdE\n>template\nA-eE\n")

        ids, seqs = read_a3m(a3m)

        @test ids == ["query", "template"]
        @test seqs == ["ACE", "A-E"]
    end

    @test AlphaConformers._auth_num("34A") == 34
    @test AlphaConformers._auth_num("A") == typemax(Int)
    @test AlphaConformers._ins_char("34A") == 'A'
    @test AlphaConformers._ins_char("34") == ' '
    @test AlphaConformers._is_peptide_resname(" ala ")
    @test AlphaConformers._is_peptide_resname("MSE")
    @test !AlphaConformers._is_peptide_resname("HOH")
    @test AlphaConformers._guess_element("") == "?"
    @test AlphaConformers._guess_element(" CA ") == "C"
    @test AlphaConformers._guess_element("ZN1") == "ZN"
    @test AlphaConformers._guess_element("XX") == "X"
end

@testitem "patch_mmcif_for_alphafold writes required template metadata" begin
    import MIToS
    import Printf

    function write_test_pdb(path)
        residues = ["ALA", "GLY", "SER"]
        serial = 1
        open(path, "w") do io
            for (resnum, resname) in enumerate(residues)
                for (atom, x, element) in
                    (("N", 0.0, "N"), ("CA", 1.0, "C"), ("C", 2.0, "C"))
                    Printf.@printf(
                        io,
                        "ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n",
                        serial,
                        atom,
                        resname,
                        resnum,
                        x + 3.0 * resnum,
                        0.0,
                        0.0,
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
        pdb = joinpath(dir, "template.pdb")
        cif = joinpath(dir, "template.cif")
        patched = joinpath(dir, "patched.cif")
        write_test_pdb(pdb)
        residues = MIToS.PDB.read_file(pdb, MIToS.PDB.PDBFile; onlyheavy = true)
        MIToS.PDB.write_file(cif, residues, MIToS.PDB.MMCIFFile)

        result = AlphaConformers.patch_mmcif_for_alphafold(
            cif,
            patched;
            entry_id = "TEST",
            force_peptide_types = true,
        )
        text = read(result, String)

        @test result == patched
        @test occursin("_entry.id", text)
        @test occursin("TEST", text)
        @test occursin("_exptl.method", text)
        @test occursin("_entity_poly_seq.entity_id", text)
        @test occursin("L-peptide linking", text)
    end
end
