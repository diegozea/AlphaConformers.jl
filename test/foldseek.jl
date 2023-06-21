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

    # Parse the output table
    df = read_foldseek_search_results(output_file)

    # Check that the results in the table are correct based on the test database 
    # tailored for this test
    @test size(df, 1) == 3 # 3 hits: 4F4J (A & B) and 1EX7
    @test df[1, :query] == "1EX6_B.pdb"
    @test df[1, :target] == "4F4J.pdb_A"
    @test df[1, :fident] > 0.9
    @test df[1, :alnlen] == 183
    @test df[1, :mismatch] == 1
    @test df[1, :gapopen] == 0
    @test df[1, :qstart] == 1
    @test df[1, :qend] == 183
    @test df[1, :tstart] == 6
    @test df[1, :tend] == 188
    @test df[1, :evalue] < 0.001
    @test df[1, :bits] == 1216

    # Clean up
    isfile(output_file) && rm(output_file)
end