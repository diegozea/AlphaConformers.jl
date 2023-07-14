@testitem "get_residues_and_sequence" begin
    query_pdb = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    query_info = get_residues_and_sequence(query_pdb)
    @test query_info.sequence == "SRPIVISGPSGTGKSTLLKKLFAEYPDSFGFSVSSTTRTPRAGEVNGKDYNF" *
        "VSVDEFKSMIKNNEFIEWAQFSGNYYGSTVASVKQVSKSGKTCILDIDMQGVKSVKAIPELNARFLFIAPPSVEDLK" * 
        "KRLEGRGTETEESINKRLSAAQAELAYAETGAHDKVIVNDDLDKAYKELKDFIFAEK"
    @test length(query_info.residues) == length(query_info.sequence)
    @test query_info.residues[1].id.name == "SER"
    @test query_info.residues[1].id.number == "201"
    @test query_info.residues[1].id.chain == "B"
    @test query_info.residues[1].id.model == "1"
    @test length(query_info.residues[1].atoms) == 6
end

@testitem "save_sequences" begin
    pdb_db = joinpath(@__DIR__, "data", "test_db")
    targets = ["4F4J.pdb", "1EX7.pdb"]
    pdb_files = [joinpath(pdb_db, target) for target in targets]
    
    mktempdir() do tmp_folder
        # The PDB files have more than one chain...
        @test_throws AssertionError save_sequences(tmp_folder, pdb_files)

        # ...so we specify the chain
        save_sequences(tmp_folder, pdb_files, chains=["A", "A"])
        filename = joinpath(tmp_folder, "sequences.fasta")
        @test isfile(filename)
        file_content = read(filename, String)
        @test occursin(">4F4J.pdb", file_content)
        @test occursin(">1EX7.pdb", file_content)
        @test length(split(file_content, "\n")) == 5 # 4 lines, but there is a \n at the end
    end
end
