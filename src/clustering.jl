# Implementation of the Hobohm algorithm I for clustering protein structures using
# the RMSD of the Cα atoms measured by USalign.

# The sorting step will not be performed, but we will use the query structure as the 
# first cluster seed.

"""
    StructuralClustering

Struct to store the results of the structural clustering. It contains the absolute paths
to the PDB files (`pdbs`), the cluster number assigned to each file (`clusters`),
the size of each cluster (`cluster_sizes`), and the total number of clusters (`nclusters`).
This struct is a subtype of `Clustering.ClusteringResult` and implements its interface by 
defining the `nclusters`, `counts`, and `assignments` methods.
"""
struct StructuralClustering <: Clustering.ClusteringResult
    pdbs::Vector{String}
    clusters::Vector{Int}
    cluster_sizes::Vector{Int}
    nclusters::Int
end

Clustering.nclusters(clustering::StructuralClustering) = clustering.nclusters
Clustering.counts(clustering::StructuralClustering) = clustering.cluster_sizes
Clustering.assignments(clustering::StructuralClustering) = clustering.clusters

"""
    _get_abspath(list_item, pdb_folder)

Helper function to get the absolute path of a file. This is useful as the `pdb_list` for 
`usalign` are relative paths to a subset of the files in `pdb_folder`. However, the first 
argument of `usalign` must be an absolute path.
"""
function _get_abspath(list_item, pdb_folder)
    if isabspath(list_item)
        list_item
    else
        joinpath(pdb_folder, list_item)
    end
end

function structural_clustering(query_pdb, pdb_folder, targets; rmsd_cutoff::Float64=1.0)
    # Ensure that the query structure is before the target structures
    pdb_list = String[query_pdb]
    append!(pdb_list, targets)
    unique!(pdb_list)
    
    # List of the absolute paths to the PDB files
    pdbs = [ _get_abspath(pdb, pdb_folder) for pdb in pdb_list ]

    N = length(pdb_list)
    clusters = zeros(Int, N)
    cluster_sizes = Int[]
    cluster_number = 0
    for i in 1:(N-1)
        if clusters[i] == 0
            cluster_number += 1
            cluster_seed = pdbs[i]
            clusters[i] = cluster_number
            push!(cluster_sizes, 1)
            
            # Only select the targets that have not been clustered yet
            pos_targets = [ pos for pos in i+1:N if clusters[pos] == 0 ]
            pdb_targets = pdb_list[pos_targets]
            # Perform the structural alignment of targets to the cluster seed
            aln_data = usalign(cluster_seed, pdb_folder, pdb_targets)
            # Find the targets that are in the cluster
            in_cluster_rows = findall(≤(rmsd_cutoff), aln_data.RMSD)
            in_cluster_pos = pos_targets[in_cluster_rows]
            # Assign the cluster number to the targets in the cluster
            clusters[in_cluster_pos] .= cluster_number
            # Update the cluster sizes
            cluster_sizes[cluster_number] += length(in_cluster_pos)
        end
    end
    if clusters[N] == 0
        cluster_number += 1
        clusters[N] = cluster_number
        push!(cluster_sizes, 1)
    end
    StructuralClustering(pdbs, clusters, cluster_sizes, cluster_number)
end