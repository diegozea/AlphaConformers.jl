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
