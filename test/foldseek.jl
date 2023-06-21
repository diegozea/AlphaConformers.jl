@testitem "Foldseek search" begin
    input_pdb = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    output_file = joinpath(@__DIR__, "data", "1EX6_B_results.m8")
    test_db_folder = joinpath(@__DIR__, "data", "test_db")

    # Throw an error if db_path is not set
    @test_throws ErrorException foldseek_search(input_pdb)

    # Throw an error if db_path is not a file
    @test_throws ErrorException foldseek_search(input_pdb, db_path=test_db_folder)

    # Run it
    foldseek_search(input_pdb, db_path=joinpath(test_db_folder, "test_db"))

    # Check that the output file exists
    @test isfile(output_file)

    # Clean up
    isfile(output_file) && rm(output_file)
end