import DataFrames
import PairwiseListMatrices
import Clustering
import Statistics

function run_tmalign_and_cluster(
    local_pdb_folder::String,
    filename_prefix::String = "tmalign",
)
    # create a temporary file for the chain list
    chain_list_path = tempname()
    open(chain_list_path, "w") do io
        for file in readdir(local_pdb_folder)
            if occursin(".pdb", file) # only select PDB files
                println(io, file)
            end
        end
    end

    # run TMalign with the -dir option:
    dir_path = abspath(local_pdb_folder)
    # add / at the end of the path if it is not there
    if dir_path[end] != '/'
        dir_path *= '/'
    end

    stdout_path = filename_prefix * ".out"
    stderr_path = filename_prefix * ".err"

    run(
        pipeline(
            `/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/TMalign -dir $dir_path $chain_list_path -fast`,
            stdout = stdout_path,
            stderr = stderr_path,
        ),
    )
    tmalign_out = read(stdout_path, String)

    tmalign_data = DataFrames.DataFrame(parse_tm_align_output(tmalign_out))

    rmsd_mat = PairwiseListMatrices.from_table(tmalign_data[:, 1:3], false)

    rmsd_clust = Clustering.hclust(rmsd_mat, linkage = :complete, branchorder = :optimal)

    rmsd_clusters = Clustering.cutree(rmsd_clust, h = 1.0)

    rmsd_chains = names(rmsd_mat)[1]

    # Return named tuple
    return (
        tmalign_data = tmalign_data,
        rmsd_mat = rmsd_mat,
        rmsd_chains = rmsd_chains,
        rmsd_clusters = rmsd_clusters,
    )
end

"""
    parse_tm_align_output(tm_align_output::String)

Parse the output of a TM-align alignment using `-dir` and return a named tuple containing
the name of Chain_1, the name of Chain_2, RMSD, Seq_ID, and the TM-scores 
(normalized by the length of Chain_1 and Chain_2) for each alignment.

# Arguments
- `tm_align_output::String`: The text output from TM-align.

# Returns
- `Vector{NamedTuple}`: A vector of named tuples, where each tuple corresponds to 
  an alignment and contains the following fields: `:chain_a`, `:chain_b`, `:RMSD`, 
  `:seq_id`, `:TM_score_a`, and `:TM_score_b`.
"""
function parse_tm_align_output(tm_align_output::String)
    alignments = split(tm_align_output, r" \*+")
    # TODO: improve this by iterating the lines and looking for Chain_1 as an alignment start
    parsed_alignments = NamedTuple{
        (:chain_a, :chain_b, :RMSD, :seq_id, :TM_score_a, :TM_score_b),
        Tuple{String,String,Vararg{Float64,4}},
    }[]

    for alignment in alignments
        if alignment != "" && occursin("Chain_1", alignment)
            chain_a = match(r"Chain_1: .*/(.*\.pdb_?[a-zA-Z0-9]?)", alignment).captures[1]
            chain_b = match(r"Chain_2: .*/(.*\.pdb_?[a-zA-Z0-9]?)", alignment).captures[1]
            RMSD = parse(Float64, match(r"RMSD=\s+(.*),", alignment).captures[1])
            seq_id = parse(
                Float64,
                match(r"Seq_ID=n_identical/n_aligned=\s+(.*)", alignment).captures[1],
            )
            TM_score_a = parse(
                Float64,
                match(r"TM-score=\s+(.*) \(if normalized by length of Chain_1", alignment).captures[1],
            )
            TM_score_b = parse(
                Float64,
                match(r"TM-score=\s+(.*) \(if normalized by length of Chain_2", alignment).captures[1],
            )
            push!(
                parsed_alignments,
                (
                    chain_a = chain_a,
                    chain_b = chain_b,
                    RMSD = RMSD,
                    seq_id = seq_id,
                    TM_score_a = TM_score_a,
                    TM_score_b = TM_score_b,
                ),
            )
        end
    end

    return parsed_alignments
end


const WORKING_DIR = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/example/clusters/cluster_3"
cd(WORKING_DIR)

"""
This function takes the path to an AlphaFold run output. 
Creates a temporal folder where the templates and models are located. 
Then run TM-align to cluster them and returns the results.
"""
function cluster_models(alphafold_output_path, templates_path)
    # Create a temporal folder
    tmp_dir = mktempdir()

    # Copy the templates and models to the temporal folder
    for file in readdir(templates_path)
        cp(joinpath(templates_path, file), joinpath(tmp_dir, file))
    end
    for file in readdir(alphafold_output_path, join = true)
        cp(file, joinpath(tmp_dir, basename(file)))
    end
    # Run TM-align and cluster the results
    return run_tmalign_and_cluster(tmp_dir, "tmalign_results")
end

# Run the function for the af2, af5 and af6 runs
templates_path = "templates"
af_results_paths = ["af2", "af5", "af6", "af7", "af8", "af9"]
af_clusters = map(af_results_paths) do path
    alphafold_output_path = joinpath(path, "predictions", "sequences", "models")
    cluster_models(alphafold_output_path, templates_path)
end

function rmsd_stats(index)
    data = af_clusters[index].tmalign_data
    sel = [
        xor(startswith(row.chain_a, "sequences"), startswith(row.chain_b, "sequences"))
        for row in eachrow(data)
    ]
    rmsd = data[sel, :RMSD]
    minimum(rmsd), Statistics.median(rmsd), maximum(rmsd)
end

# julia> rmsd_stats(1) # af2
# (2.88, 4.685, 6.41)

# julia> rmsd_stats(2) # af5
# (2.65, 4.315, 6.49)

# julia> rmsd_stats(3) # af6
# (2.66, 4.18, 6.28)

# julia> rmsd_stats(4) # af7
# (2.12, 4.475, 6.7)

# julia> rmsd_stats(5) # af8  <--- best (only models using templates)
# (2.53, 4.675, 6.57)

# julia> rmsd_stats(6) # af9
# (2.76, 4.55, 6.42)


# 1217  /opt/alphafold/runcolabfold.py sequences.a3m af2 --msa-input --custom-template-path templates/
# 1232  /opt/alphafold/runcolabfold.py sequences.a3m af5 --use-templates 1 --msa-input --custom-template-path templates/
# 1233  /opt/alphafold/runcolabfold.py sequences.a3m af6 --use-templates 1 --msa-input
# 1249  /opt/alphafold/runcolabfold.py sequences.a3m af7 --use-templates 1 --msa-input --custom-template-path templates/ --num-seeds 5 --use-dropout
# 1258  /opt/alphafold/runcolabfold.py sequences.a3m af8 --use-templates 1 --msa-input --custom-template-path templates/ --num-seeds 5 --use-dropout --num-models 2
# 1270  /opt/alphafold/runcolabfold.py sequences.a3m af9 --use-templates 1 --msa-input --custom-template-path templates/ --num-seeds 5 --use-dropout --num-models 2 --num-ensemble 5

# I have run af8 for the other 3 clusters
# So, let's check the RMSD for the models

cd("/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/example/clusters/")
for folder in ["cluster_$i" for i in [1, 2, 3, 4]]
    cd(folder)
    templates_path = "templates"
    alphafold_output_path = "af8/predictions/sequences/models"
    af8_cluster = cluster_models(alphafold_output_path, templates_path)
    data = af8_cluster.tmalign_data
    sel = [
        xor(startswith(row.chain_a, "sequences"), startswith(row.chain_b, "sequences"))
        for row in eachrow(data)
    ]
    rmsd = data[sel, :RMSD]
    println("Cluster $folder")
    println(
        "min: $(minimum(rmsd)), median: $(Statistics.median(rmsd)), max: $(maximum(rmsd))",
    )
    cd("..")
end

# Cluster cluster_1
# min: 1.29, median: 3.58, max: 6.66
# Cluster cluster_2
# min: 2.61, median: 4.885, max: 6.7
# Cluster cluster_3
# min: 2.53, median: 4.675, max: 6.57
# Cluster cluster_4
# min: 4.76, median: 5.795, max: 6.39
