# Implementation of the Hobohm algorithm I for clustering protein structures using
# the RMSD of the Cα atoms measured by USalign.

# The sorting step will not be performed, but we will use the query structure as the 
# first cluster seed.


function structure_clustering(query_pdb, pdb_folder, targets)
    pdb_list = String[query_pdb]
    append!(pdb_list, targets)
    unique!(pdb_list)
    N = length(pdb_list)
    clusters = zeros(Int, )
    cluster_number = 0
    for i in 1:(N-1)
        if clusters[i] == 0
            cluster_number += 1
            query = pdb[i]
            targets = pdb[i+1:N]
            # clustering ... to continue
        end
    end
    if clusters[N] == 0
        cluster_number += 1
        clusters[N] = cluster_number
    end
    (pdb=pdb_list, cluster=clusters)
end