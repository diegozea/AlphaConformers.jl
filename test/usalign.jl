@testitem "USalign" begin
    import DataFrames

    query_pdb = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    pdb_db = joinpath(@__DIR__, "data", "test_db")
    targets = ["4F4J.pdb_A", "1EX7.pdb"]
    
    mktempdir() do tmp_folder
        pdb_folder = joinpath(tmp_folder, "test_folder")
        create_pdb_folder(targets, pdb_folder, pdb_db=pdb_db)
        # Check usalign (pairs)
        for target in targets
            target_path = joinpath(pdb_folder, target)
            one2one = usalign(query_pdb, target_path)
            # Check that the first column has been renamed
            # Check that the DataFrame has only one row
            @test startswith(only(one2one.PDBchain1), query_pdb) # query_pdb:chain
            @test startswith(only(one2one.PDBchain2), target_path) # target_pdb:chain
        end
        # Check usalign (one to many)
        one2many = usalign(query_pdb, pdb_folder, targets)
        @test DataFrames.nrow(one2many) == 2
        @test startswith(only(unique(one2many.PDBchain1)), query_pdb)
        # Check the @warn message when the file do not exist (_save_list)
        @test_logs (:warn, r"it won't be included in the comparison") usalign(query_pdb, pdb_folder, ["1EX7.pdb" "1PBE.pdb"])

        # Check usalign keyword arguments
        target_pdb = joinpath(pdb_folder, "4F4J.pdb_A")
        one2one_default = usalign(query_pdb, target_pdb)
        one2one_slow = usalign(query_pdb, target_pdb, TMcut=nothing, fast=false, 
            extra_parameters=`-ter 1`) # -ter has no effect as both files have one chain
        @test one2one_slow == one2one_default
        # When USalign fail to give the correct output, reading the DataFrame fails
        @test_throws ArgumentError usalign(query_pdb, target_pdb, outfmt=0)
        # Check that extra_parameters is not ignored by giving a wrong parameter
        @test_throws ProcessFailedException usalign(query_pdb, target_pdb, extra_parameters=`-dir2`)
    end
end