@testitem "UniProt PDB mapping" begin
    import DataFrames

    mktempdir() do tmp_folder
        data_path = joinpath(tmp_folder, "pdb_chain_uniprot.csv.gz")
        @test !isfile(data_path)
        data = get_uniprot_mapping(tmp_folder)
        @test isfile(data_path)

        # Check the DataFrame content
        @test isa(data, DataFrames.DataFrame)
        for col_name in ["PDB", "CHAIN", "SP_PRIMARY"]
            @test col_name in names(data)
        end
        @test "101m" in data.PDB
        @test "A" in data.CHAIN
        @test "P02185" in data.SP_PRIMARY
        @test !("None" in skipmissing(data.PDB_END))
        @test sum(ismissing.(data.PDB_END)) != 0

        # Check the function get_uniprot_acc
        @test get_uniprot_acc(data, "101m") == ["P02185"]

        # Check the function get_pdb_codes
        pdbs = get_pdb_codes(data, "P02185") # == [("101M", "A")]
        @test isa(pdbs, DataFrames.DataFrame)
        @test DataFrames.ncol(pdbs) == 2
        @test names(pdbs) == ["PDB", "CHAIN"]
        @test "101M" in pdbs.PDB
        @test "A" in pdbs.CHAIN
    end
end

@testitem "Updating FoldSeek search results" begin
    import DataFrames

    # Run FoldSeek and download the SIFTS file to generate the test data
    input_pdb = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    output_file = joinpath(@__DIR__, "data", "1EX6_B_results.m8")
    test_db_folder = joinpath(@__DIR__, "data", "test_db")
    foldseek_search(input_pdb, db_path=joinpath(test_db_folder, "test_db"))
    search_results = read_foldseek_search_results(output_file)
    sifts_uniprot_mapping = get_uniprot_mapping()

    # Test the list_known_conformations function
    pdb_chains = list_known_conformations(search_results, sifts_uniprot_mapping)
    for item in pdb_chains
        @test !in(item, search_results.target)
    end
    for pdb_chain in ["1EX6.pdb_B", "1EX6.pdb_A", "1GKY.pdb_A"]
        @test in(pdb_chain, pdb_chains)
    end

    # Test the delete_query_from_target! function
    
    # 1PBE A is another protein (P00438)
    @test DataFrames.nrow(delete_query_from_target!(deepcopy(search_results),  
        sifts_uniprot_mapping, "1PBE", "A")) == 3
    # 4F4J and 1EX7 are known conformations of 1EX6 B (P15454)
    @test isempty(delete_query_from_target!(deepcopy(search_results), 
        sifts_uniprot_mapping, "1EX6", "B"))
end