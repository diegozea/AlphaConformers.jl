@testitem "Foldseek search" begin
    input_pdb = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    test_db_folder = joinpath(@__DIR__, "data", "test_db")

    # Throw an error if db_path is not set
    @test_throws ErrorException foldseek_search(input_pdb)

    # Throw an error if db_path is not a file
    @test_throws ErrorException foldseek_search(input_pdb, db_path=test_db_folder)

    # Run it
    output = foldseek_search(input_pdb, db_path=joinpath(test_db_folder, "test_db"))

    # Check that the output file exists
    @test output == joinpath(@__DIR__, "data", "1EX6_B_results.m8")
    @test isfile(output)

    # Parse the output table
    df = read_foldseek_search_results(output)

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
    isfile(output) && rm(output)
end

@testitem "Foldseek aln structures" begin
    input_pdb = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    test_db_folder = joinpath(@__DIR__, "data", "test_db")

    # Run it
    output = foldseek_search(input_pdb, db_path=joinpath(test_db_folder, "test_db"), format_mode=5)

    # Check that the output folder exists
    @test output == joinpath(@__DIR__, "data", "aligned_structures")
    @test isdir(output)

    # Check that the output folder contains the expected files
    files = readdir(output)
    @test length(files) == 3
    @test all(startswith("aln_"), files)

    # Clean up
    isdir(output) && rm(output; recursive=true)
end

@testitem "run_foldseek" begin
    input_pdb = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    test_db = joinpath(@__DIR__, "data", "test_db", "test_db")

    # Run it
    output = run_foldseek(input_pdb, db_path=test_db)
    
    println(output)
end
