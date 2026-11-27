using TestItems

@testitem "Clustering helpers" begin
    import OrderedCollections

    aln = [('A', 'A'), ('B', '-'), ('C', 'D'), ('-', 'E')]
    coverage, identity = AlphaConformers.coverage_and_identity(aln)
    @test coverage == 2 / 3
    @test identity == 1 / 2

    labels =
        OrderedCollections.OrderedDict("target_a" => 1, "target_b" => 1, "target_c" => 2)
    cluster2targets, nclusters = AlphaConformers.get_cluster2targets(labels, [1, 2])
    @test nclusters == 2
    @test cluster2targets[1] == ["target_a", "target_b"]
    @test cluster2targets[2] == ["target_c"]

    target2sequence =
        Dict("target_a" => "seq_a", "target_b" => "seq_b", "target_c" => "seq_c")
    cluster2seqnames =
        AlphaConformers.get_cluster2seqnames(cluster2targets, target2sequence)
    @test cluster2seqnames[1] == ["seq_a", "seq_b"]
    @test cluster2seqnames[2] == ["seq_c"]

    duplicate_targets = [("target_a", 1), ("target_a", 1), ("target_b", 11)]
    deduped, ndeduped = AlphaConformers.get_cluster2targets(duplicate_targets, [1, 11])
    @test ndeduped == 2
    @test deduped[1] == ["target_a"]

    concatenated, nconcatenated =
        AlphaConformers.get_cluster2targets_concatenate(duplicate_targets, [1, 11])
    @test nconcatenated == 2
    @test concatenated[1] == ["target_a"]
    @test concatenated[2] == ["target_b"]

    @test AlphaConformers._seq_name_to_key("target\t12\t0.5") == ("target", 12, 0.5)
    @test AlphaConformers._find_duplicates(["a", "b", "a"]) == Set(["a"])

    missing_sequence = AlphaConformers.get_cluster2seqnames(
        OrderedCollections.OrderedDict(1 => ["target_a", "missing"]),
        Dict("target_a" => "seq_a"),
    )
    @test missing_sequence[1] == ["seq_a"]
end

@testitem "clustering folder helpers handle existing folders safely" begin
    import MIToS
    import OrderedCollections

    mktempdir() do dir
        out_folder = joinpath(dir, "nested", "run")
        c1 = joinpath(out_folder, "cluster_1")
        mkpath(joinpath(c1, "templates_complete"))

        msa = MIToS.MSA.parse_file(IOBuffer(">1ABC_A\nACD\n>2DEF_A\nA-D\n"), MIToS.MSA.FASTA)
        cl2msa = OrderedCollections.OrderedDict(1 => msa)
        cl2pdb = OrderedCollections.OrderedDict(
            1 => Dict{String,Vector{MIToS.PDB.PDBResidue}}(),
        )

        @test AlphaConformers.create_folder_structure_hobohm(
            1,
            cl2msa,
            cl2pdb;
            out_folder = out_folder,
        ) == out_folder
        @test isfile(joinpath(c1, "sequences.a3m"))
        @test isfile(joinpath(out_folder, "calpha_template.csv"))

        # download_template and clean_msa_template_names reuse/skip existing folders cleanly
        AlphaConformers.download_template(c1, String[])
        write(joinpath(c1, "templates_complete", "2def.cif"), "data_2def\n")
        AlphaConformers.clean_msa_template_names(1, out_folder)
        @test isfile(joinpath(c1, "templates_adaptative", "2def.cif"))
        @test isfile(joinpath(c1, "corresponding_name_template.csv"))
        redirect_stdout(devnull) do
            AlphaConformers.clean_msa_template_names(1, out_folder)
        end
        @test isfile(joinpath(c1, "templates_adaptative", "2def.cif"))
    end
end
