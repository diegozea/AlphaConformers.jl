import Pkg
# Pkg.activate(abspath(@__DIR__, ".."))
Pkg.activate("/home/diego.zea/.julia/dev/AlphaConformers/")

import DataFrames
import Statistics
import Plots
import StatsPlots
import PairwiseListMatrices
import Clustering
import MIToS

Plots.plotly()

const WORKING_DIR = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/from_apo"

const folders = ["1AEL_A", "1BV2_A", "1C54_A", "1EAL_A"]
# 1AEL: NMR
# 1BV2: NMR
# 1C54: NMR
# 1EAL: NMR 
# ... all the runned examples have the bug with the NMR sequences.

function progress(folder)
    results = NamedTuple{(:folder, :subfolder, :finished),Tuple{String,String,Bool}}[]
    cd(WORKING_DIR)
    cd(folder)
    cd("clusters")
    subfolders = filter(f -> isdir(f) && startswith(f, "cluster_"), readdir())
    for subfolder in subfolders
        cd(subfolder)
        progress = readlines("af/progress.log")
        finished = occursin("Job finished", progress[end])
        push!(results, (; folder, subfolder, finished))
        cd("..")
    end
    cd("..")
    results
end

results = map(folders) do folder
    progress(folder)
end

data = DataFrames.DataFrame(vcat(results...))

finished =
    DataFrames.combine(DataFrames.groupby(data, "folder"), "finished" => Statistics.mean)

const finished_folders = finished.folder[finished.finished_mean .== 1.0]

# ------------------------------------------------------------------------------------------
# From compare_runs.jl


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

function rmsd_stats(index)
    data = af_clusters[index].tmalign_data
    sel = [
        xor(startswith(row.chain_a, "sequences"), startswith(row.chain_b, "sequences"))
        for row in eachrow(data)
    ]
    rmsd = data[sel, :RMSD]
    minimum(rmsd), Statistics.median(rmsd), maximum(rmsd)
end

# ------------------------------------------------------------------------------------------

function analyze_clusters(WORKING_DIR::String, folder::String)
    results = []
    for subfolder in readdir(joinpath(WORKING_DIR, folder, "clusters"), join = true)
        cd(subfolder)
        templates_path = "templates"
        alphafold_output_path = "af/predictions/sequences/models"
        af_cluster = cluster_models(alphafold_output_path, templates_path)
        push!(results, (; folder, af_cluster))
        data = af_cluster.tmalign_data
        sel = [
            xor(startswith(row.chain_a, "sequences"), startswith(row.chain_b, "sequences")) for row in eachrow(data)
        ]
        rmsd = data[sel, :RMSD]
        println("Cluster $subfolder")
        println(
            "min: $(minimum(rmsd)), median: $(Statistics.median(rmsd)), max: $(maximum(rmsd))",
        )
        cd("..")
    end
    results
end

#=
a = analyze_clusters(WORKING_DIR, finished_folders[1])
b = analyze_clusters(WORKING_DIR, finished_folders[2])
=#

using Serialization

function save_objects_to_disk(objects::Dict{String,Any})
    for (filename, obj) in objects
        open(filename, "w") do f
            serialize(f, obj)
        end
    end
end

cd(WORKING_DIR)

#save_objects_to_disk(Dict{String, Any}("a.bin" => a, "b.bin" => b))

function rmsd_stats_cluster(af_cluster)
    data = af_cluster.tmalign_data
    sel = [
        xor(startswith(row.chain_a, "sequences"), startswith(row.chain_b, "sequences"))
        for row in eachrow(data)
    ]
    rmsd = data[sel, :RMSD]
    println(
        "min: $(minimum(rmsd)), median: $(Statistics.median(rmsd)), max: $(maximum(rmsd))",
    )
end

#=
for elem in a
    println(elem.folder)
    rmsd_stats_cluster(elem.af_cluster)
end

for elem in b
    println(elem.folder)
    rmsd_stats_cluster(elem.af_cluster)
end
=#

using Serialization

function load_object_from_disk(filename::String)
    open(filename, "r") do f
        return deserialize(f)
    end
end

cd(WORKING_DIR)
a = load_object_from_disk("a.bin")
b = load_object_from_disk("b.bin")

# 1BV2-3_A;1UVB_A;1.9621;1.3671;1.4394;96.94;91;6170;1.13;Holo;88.97;3;domain movements;4.95
# 1C54-19_A;1RGH_B;1.4925;1.515;0.4235;97.29;96;5679;0.42;Holo;97.29;1;loop movements;4.79

#=
# Set up the known conformations
MIToS.PDB.downloadpdb("1BV2", filename="1BV2.pdb.gz", format=MIToS.PDB.PDBFile)
pdb_a1 = read("1BV2.pdb.gz", MIToS.PDB.PDBFile, model="1", chain="A", onlyheavy=true, occupancyfilter=true)
write("1BV2.pdb", pdb_a1, MIToS.PDB.PDBFile)
MIToS.PDB.downloadpdb("1UVB", filename="1UVB.pdb.gz", format=MIToS.PDB.PDBFile)
pdb_a2 = read("1UVB.pdb.gz", MIToS.PDB.PDBFile, model="1", chain="A", onlyheavy=true, occupancyfilter=true)
write("1UVB.pdb", pdb_a2, MIToS.PDB.PDBFile)
MIToS.PDB.downloadpdb("1C54", filename="1C54.pdb.gz", format=MIToS.PDB.PDBFile)
pdb_b1 = read("1C54.pdb.gz", MIToS.PDB.PDBFile, model="1", chain="A", onlyheavy=true, occupancyfilter=true)
write("1C54.pdb", pdb_b1, MIToS.PDB.PDBFile)
MIToS.PDB.downloadpdb("1RGH", filename="1RGH.pdb.gz", format=MIToS.PDB.PDBFile)
pdb_b2 = read("1RGH.pdb.gz", MIToS.PDB.PDBFile, model="1", chain="B", onlyheavy=true, occupancyfilter=true)
write("1RGH.pdb", pdb_b2, MIToS.PDB.PDBFile)
=#

function get_models(folder)
    clusters = readdir(joinpath(WORKING_DIR, folder, "clusters"), join = true)
    list_models = map(clusters) do cluster
        models = readdir(joinpath(cluster, "af/predictions/sequences/models"), join = true)
    end
    vcat(list_models...)
end

function run_tmalign_and_cluster_models(folder, pdb_a, pdb_b)
    models = get_models(folder)
    tmp_dir = mktempdir()
    cp(pdb_a, joinpath(tmp_dir, basename(pdb_a)))
    cp(pdb_b, joinpath(tmp_dir, basename(pdb_b)))
    for model in models
        cp(model, joinpath(tmp_dir, basename(model)))
    end
    return run_tmalign_and_cluster(tmp_dir, "tmalign_$folder")
end

#=
data_a = run_tmalign_and_cluster_models("1BV2_A", "1BV2.pdb", "1UVB.pdb")
data_b = run_tmalign_and_cluster_models("1C54_A", "1C54.pdb", "1RGH.pdb")
save_objects_to_disk(Dict{String, Any}("data_a.bin" => data_a, "data_b.bin" => data_b))
=#

cd(WORKING_DIR)
data_a = load_object_from_disk("data_a.bin")
data_b = load_object_from_disk("data_b.bin")

function plot_superimposition(pdb_file_a, pdb_file_b)
    # Read the pdb files
    res_a = MIToS.PDB.read(pdb_file_a, MIToS.PDB.PDBFile, atomname = "CA")
    res_b = MIToS.PDB.read(pdb_file_b, MIToS.PDB.PDBFile, atomname = "CA")
    @show length(res_a)
    @show length(res_b)
    # Superimpose the chains
    superimposed_A, superimposed_B, RMSD = MIToS.PDB.superimpose(res_a, res_b)

    # Print the RMSD value
    println("RMSD: ", RMSD)

    # Plot the superimposed chains
    Plots.scatter3d(superimposed_A, label = "A", alpha = 0.5)
    Plots.scatter3d!(superimposed_B, label = "B", alpha = 0.5)
end

function find_file_path(filename::AbstractString, search_directory::AbstractString)
    command = `find $search_directory -name $filename`
    strip(read(command, String))
end

# Usage example
filename = "sequences_unrelaxed_rank_006_alphafold2_ptm_model_1_seed_41003.pdb"
search_directory = pwd()
file_path = find_file_path(filename, search_directory)
