@testitem "Clustering" begin
    import Clustering

    query_pdb = joinpath(@__DIR__, "data", "1EX6_B.pdb")
    pdb_db = joinpath(@__DIR__, "data", "test_db")
    targets = ["4F4J.pdb", "1EX7.pdb"]
    
    mktempdir() do tmp_folder
        pdb_folder = joinpath(tmp_folder, "test_folder")
        create_pdb_folder(targets, pdb_folder, pdb_db=pdb_db)
        
        # Each structure is in a different cluster
        clusters = structural_clustering(query_pdb, pdb_folder, targets) # rmsd_cutoff=1.0
        @test Clustering.nclusters(clusters) == 3
        @test Clustering.counts(clusters) == [1, 1, 1]
        @test Clustering.assignments(clusters) == [1, 2, 3]
        @test [ basename.(l) for l in get_clustered_pdbs(clusters) ] == [
            ["1EX6_B.pdb"], ["4F4J.pdb"], ["1EX7.pdb"]]

        # The structures are clustered in two groups
        clusters_3 = structural_clustering(query_pdb, pdb_folder, targets, rmsd_cutoff=3.0)
        @test Clustering.nclusters(clusters_3) == 2
        @test Clustering.counts(clusters_3) == [2, 1]
        @test Clustering.assignments(clusters_3) == [1, 1, 2]
        @test [ basename.(l) for l in get_clustered_pdbs(clusters_3) ] == [
            ["1EX6_B.pdb", "4F4J.pdb"], ["1EX7.pdb"]]

        # All structures are in the same cluster
        clusters_100 = structural_clustering(query_pdb, pdb_folder, targets, 
            rmsd_cutoff=100.0)
        @test Clustering.nclusters(clusters_100) == 1
        @test Clustering.counts(clusters_100) == [3]
        @test Clustering.assignments(clusters_100) == [1, 1, 1]
        @test [ basename.(l) for l in get_clustered_pdbs(clusters_100) ] == [
            ["1EX6_B.pdb", "4F4J.pdb", "1EX7.pdb"]]

        for cluster_result in [clusters, clusters_3, clusters_100]
            # Check that the query structure is the first one and that the targets are
            # in the same order as in the input
            @test basename.(cluster_result.pdbs) == ["1EX6_B.pdb", "4F4J.pdb", "1EX7.pdb"]
            # Check that the vectors have the correct length
            @test length(cluster_result.pdbs) == length(cluster_result.clusters)
            @test length(cluster_result.cluster_sizes) == cluster_result.nclusters
        end
    end
end