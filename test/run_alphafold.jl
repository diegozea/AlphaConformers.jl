@testitem "run_alphafold" begin
    
    mktempdir() do temp_folder
        # test the error when the colabfold_path is not defined
        @test_throws ErrorException run_alphafold(temp_folder)
        # test the error when no cluster_* folders are found
        @test_throws ErrorException run_alphafold(temp_folder, colabfold_path="echo")
    end
end