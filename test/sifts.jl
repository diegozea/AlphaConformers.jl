@testitem "PDB folders" begin
    import DataFrames

    mktempdir() do tmp_folder
        data_path = joinpath(tmp_folder, "pdb_chain_uniprot.csv.gz")
        @test !isfile(data_path)
        data = get_uniprot_mapping(tmp_folder)
        @test isfile(data_path)
        @test isa(data, DataFrames.DataFrame)
        for col_name in ["PDB", "CHAIN", "SP_PRIMARY"]
            @test col_name in names(data)
        end
        @test "101m" in data.PDB
        @test "A" in data.CHAIN
        @test "P02185" in data.SP_PRIMARY
        @test !("None" in skipmissing(data.PDB_END))
        @test sum(ismissing.(data.PDB_END)) != 0
    end
end