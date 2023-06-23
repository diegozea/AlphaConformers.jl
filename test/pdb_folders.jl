@testitem "PDB folders" begin
    import DataFrames
    import MIToS

    mktempdir() do tmp_folder
        folder_path = joinpath(tmp_folder, "test_folder")
   
        targets = Set(["1EX7.pdb", "4F4J.pdb_B", # Examples from a local database
            "1lvg", # 1 chain, without extension
            "6nui_A", # Example from the webserver; 6nui has 1 chains and 20 models (NMR)
            "6mfu_C", # Example from the webserver; 6mfu has 4 chains
            ])

        # Check if a folder is created when the function is called with a Set
        created_folder_path = create_pdb_folder(targets, folder_path)
        @test abspath(created_folder_path) == abspath(folder_path)
        @test isdir(created_folder_path)

        # Check if the correct number of files is created
        @test length(readdir(created_folder_path)) == length(targets)

        # Check that "1lvg" is in the folder
        @test isfile(joinpath(created_folder_path, "1lvg"))

        # Check if the correct files are created
        for target in targets
            file_path = joinpath(created_folder_path, target)
            @test isfile(file_path)
            res = read(file_path, MIToS.PDB.PDBFile)
            # Check that there is only one chain
            @test length(unique([ r.id.chain for r in res])) == 1
        end

        # Check if the function works with a DataFrame
        df = DataFrames.DataFrame(target = ["1EX7.pdb", "4F4J.pdb_B"])
        created_folder_path = create_pdb_folder(df, folder_path)
        @test isdir(created_folder_path)
        @test abspath(created_folder_path) == abspath(folder_path)
        @test length(readdir(created_folder_path)) == DataFrames.nrow(df)

        # Check that "1lvg" is not in the folder as the previous folder has been deleted
        @test !isfile(joinpath(created_folder_path, "1lvg"))

        # Check if the function works with an array
        arr = ["1EX7.pdb", "4F4J.pdb_B"]
        created_folder_path = create_pdb_folder(arr, folder_path)
        @test isdir(created_folder_path)
        @test abspath(created_folder_path) == abspath(folder_path)
        @test length(readdir(created_folder_path)) == length(arr)
    end
end
