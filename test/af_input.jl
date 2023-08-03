@testitem "_read_pdb_chain" begin
    pdb_file = joinpath(@__DIR__, "data", "7ADD.pdb.gz")

    # auth chain is lowercase in 7ADD; test that it is read correctly
    lowercase_chain = AlphaConformers._read_pdb_chain(pdb_file, "e")
    uppercase_chain = AlphaConformers._read_pdb_chain(pdb_file, "E")
    @test lowercase_chain == uppercase_chain
end

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

@testitem "align and clean sequences" begin
    import MIToS

    pdb_db = joinpath(@__DIR__, "data", "test_db")
    targets = ["4F4J.pdb", "1EX7.pdb"]
    pdb_files = [joinpath(pdb_db, target) for target in targets]

    # Test the align_sequences function
    @test_throws AssertionError align_sequences("non_existent_file.fasta")

    msa = mktempdir() do tmp_folder
        save_sequences(tmp_folder, pdb_files, chains=["A", "A"])
        filename = joinpath(tmp_folder, "sequences.fasta")
        align_sequences(filename)
    end
    
    @test MIToS.MSA.nsequences(msa) == 2
    @test MIToS.MSA.sequencenames(msa) == ["4F4J.pdb", "1EX7.pdb"] # keep the original order
    @test MIToS.MSA.stringsequence(msa, 1) == "FQGSMSRPIVISGPSGTGKSTLLKKLFAEYPDSFG" * 
        "FSVPSTTRTPRAGEVNGKDYNFVSVDEFKSMIKNNEFIEWAQFSGNYYGSTVASVKQVSKSGKTCILDIDMQG" * 
        "VKSVKAIPELNARFLFIAPPSVEDLKKRLEGRGTETEESINKRLSAAQAELAYAETGAHDKVIVNDDLDKAYK" *
        "ELKDFIFA--"
    
    # Test the clean_msa function
    cleaned_msa = clean_msa(msa)
    @test MIToS.MSA.nsequences(cleaned_msa) == 2 # keep both sequences
    @test MIToS.MSA.sequencenames(cleaned_msa) == ["4F4J.pdb", "1EX7.pdb"]
    @test MIToS.MSA.stringsequence(cleaned_msa, 1) == "FQGSMSRPIVISGPSGTGKSTLLKKLF" * 
        "AEYPDSFGFSVPSTTRTPRAGEVNGKDYNFVSVDEFKSMIKNNEFIEWAQFSGNYYGSTVASVKQVSKSGKTC" * 
        "ILDIDMQGVKSVKAIPELNARFLFIAPPSVEDLKKRLEGRGTETEESINKRLSAAQAELAYAETGAHDKVIVN" *
        "DDLDKAYKELKDFIFA" # no gaps in the first sequence
end

@testitem "create_pdb_lists" begin
    query = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    pdb_db = joinpath(@__DIR__, "data", "test_db")
    targets = ["4F4J.pdb_A", "1EX7.pdb"]
    pdb_files = [joinpath(pdb_db, target) for target in targets]

    ref_first = pushfirst!(deepcopy(pdb_files), query)

    # Do nothing if the query is the first pdb file on the list
    ref_first_paths = create_pdb_lists(query, "B", "1", ref_first, ["B", "A", "A"], ["1", "1", "1"])
    @test ref_first_paths.pdb_files == abspath.(ref_first)
    @test ref_first_paths.chains == ["B", "A", "A"]
    @test ref_first_paths.models == ["1", "1", "1"]

    # Throw an error if the query is on the pdb list, but not in chains and models
    @test_throws AssertionError create_pdb_lists(query, "B", "1", ref_first, ["A", "A"], ["1", "1"])

    # Add the reference if it is not on the list
    paths = create_pdb_lists(query, "B", "1", pdb_files, ["A", "A"], ["1", "1"])
    @test paths.pdb_files == abspath.(ref_first)
    @test paths.chains == ["B", "A", "A"]
    @test paths.models == ["1", "1", "1"]

    # Move the reference to the first position if it is on the list
    ref_last = push!(deepcopy(pdb_files), query)

    ref_last_paths = create_pdb_lists(query, "B", "1", ref_last, ["A", "A", "B"], ["1", "1", "1"])
    @test ref_last_paths.pdb_files == abspath.(ref_first[[1, 3, 2]]) # 2 was the first
    @test ref_last_paths.chains == ["B", "A", "A"]
    @test ref_last_paths.models == ["1", "1", "1"]
end

@testitem "create_msa_and_templates" begin
    import MIToS

    query = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    pdb_db = joinpath(@__DIR__, "data", "test_db")
    targets = ["4F4J.pdb", "1EX7.pdb"]
    pdb_files = [joinpath(pdb_db, target) for target in targets]

    ref_pdb = query
    ref_chain = "B"
    ref_model = "1"
    chains = ["A", "A"]
    models = ["1", "1"]

    mktempdir() do cluster_folder
        output = create_msa_and_templates(cluster_folder, ref_pdb, ref_chain, ref_model, pdb_files, chains, models)
        # test the returned values
        @test MIToS.MSA.sequencenames(output.msa) == ["1EX6_B.pdb", "4F4J.pdb", "1EX7.pdb"]
        @test output.pdb_files[1] == abspath(query)
        @test length(output.pdb_files) == 3
        @test output.chains == ["B", "A", "A"]
        @test output.models == ["1", "1", "1"]
        
        # test the created files
        @test isfile(joinpath(cluster_folder, "sequences.fasta"))
        @test isfile(joinpath(cluster_folder, "sequences.a3m"))
        @test isdir(joinpath(cluster_folder, "templates"))
        @test length(readdir(joinpath(cluster_folder, "templates"))) == 2
    end
end

@testitem "create_alpha_fold_inputs" begin
    query = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    test_db_folder = joinpath(@__DIR__, "data", "test_db", "test_db")

    # testing = false
    mktempdir() do output_folder
        paths = create_alpha_fold_inputs(output_folder, query, "B", "1",
            foldseek_db=test_db_folder)
        @test isdir(paths.clusters)
        @test isdir(paths.pdb)
        cluster_folders = filter!(
            dir -> occursin("cluster_", dir), 
            readdir(paths.clusters, join=true))
        @test !isempty(cluster_folders)
        for folder in cluster_folders
            @test isfile(joinpath(folder, "sequences.fasta"))
            @test isfile(joinpath(folder, "sequences.a3m"))
            template_folder = joinpath(folder, "templates")
            @test isdir(template_folder)
            @test !isempty(readdir(template_folder))
        end
        @test "1EX6.pdb_A" in readdir(paths.pdb)
        @test "1EX6.pdb_B" in readdir(paths.pdb)
        @test "1EX7.pdb" in readdir(paths.pdb)
        @test "4F4J.pdb_A" in readdir(paths.pdb)
    end

    # testing = true
    mktempdir() do output_folder
        # 4F4J and 1EX7 are known conformations of 1EX6
        @test_throws ErrorException create_alpha_fold_inputs(output_folder, query, 
            "B", "1", foldseek_db=test_db_folder, testing=true)
        # Error is thrown before creating the pdb folder
        @test !isdir(joinpath(output_folder, "clusters", "pdb"))
        # But after creating the clusters folder
        @test isdir(joinpath(output_folder, "clusters"))
    end
end
