#!/stockage/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#SBATCH --nodelist=node48
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2
#SBATCH --output=output_AlphaConformers.jl.o%j.out
=#

import Pkg
Pkg.activate("/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/Clean_AlphaConformers/scripts/update_MIToS_321")


# Analyse output of AlphaConformer
using MIToS
using MIToS.MSA
using MIToS.PDB
using DataFrames, CSV
using Plots
using AlphaConformers
using Glob
using Clustering: randindex
using Glob
using Statistics
using ArgParse

################################################# FUNCTIONS #############################################

"""
Compare all the model predicted by AlphaFold to the apo and holo shape
Use the function structural_alignment from AlphaConformers
"""
#Calculate the RMSD 
function compare_model(model_files::Vector,holo_pdb::String,apo_pdb::String,folder_path::String,folder_model::String)
    results = DataFrame(Model=String[], RMSD_Holo=Float64[], RMSD_Apo=Float64[])
    for model in model_files
        println(model)

        structure_model = read_file(model, PDBFile,group="ATOM")
        structure_apo = read_file(folder_path*apo_pdb, PDBFile,group="ATOM")
        structure_holo = read_file(folder_path*holo_pdb, PDBFile,group="ATOM")

        _,_,_,rmsd_holo,coverage,_=AlphaConformers.structural_alignment(structure_model, structure_holo)
        _,_,_,rmsd_apo,coverage,_ =AlphaConformers.structural_alignment(structure_model, structure_apo)
        println(basename(model))
        push!(results, (basename(model), rmsd_holo, rmsd_apo))
        
    end
    return results 
end

"""
Get all the model predicted by ALphaFond after AlphaConformers pipeline 
Need the path to the output folder from ALphaConformer 
Output the path from all the model 

Further could use pLDDT to color the plot by confidence 
"""
#Get all thne model predicted by AlphaFold
function find_model(folder_model::String)
    model_files = String[]
    # Found the folder that start by "cluster"
    cluster_dirs = filter(isdir, glob("cluster*", folder_model))
    for cluster in cluster_dirs
        models_path = joinpath(cluster, "af/predictions/sequences/models")  #  Path for the model
        scores_file = joinpath(cluster, "af/predictions/scores.tsv")  #Path for the pLDDT
        println(models_path)
        if isdir(models_path)   # Check if the path is correct
            df_scores = CSV.read(scores_file, DataFrame, delim='\t')

            max_pLDDT = maximum(df_scores.pLDDT)
            #best_models = df_scores[df_scores.pLDDT .== max_pLDDT, :]
            #println(best_models)
            pdb_files = glob("sequences_unrelaxed_rank_*_alphafold2_ptm_model_*_seed_*.pdb", models_path)
            println(length(pdb_files))
            append!(model_files, pdb_files)
        else 
            println("Path not found")
        end
    end
    return model_files
end

"""
    analyse_output(query, folder_path, folder_model, objectif)

Compare the models predicted by AlphaFold to the apo and holo conformations,
then save a CSV and a scatter plot in the AlphaConformers output folder.

Uses `find_model` to locate all predicted PDB models and `compare_model` to
compute RMSD against both the apo and holo reference structures.

# Arguments
- `query::String`        : Name of the PDB entry being analysed (e.g. `"6akl"`).
- `folder_path::String`  : Path to the folder containing the apo and holo reference files.
- `folder_model::String` : Path to the AlphaConformers output folder.
- `objectif::String`     : Basename used to locate `<objectif>_apo.pdb` and `<objectif>_holo.pdb`.

# Outputs
- A CSV file : `<folder_model>/rmsd_results_<query>.csv`
- A PNG plot : `<folder_model>/rmsd_scatter_<query>.png`

# Notes
- TODO: color scatter points by pLDDT confidence score.
- To compare two AlphaConformers runs (number of clusters, MSA lines, templates),
  call this function on each folder separately and compare the resulting CSVs.

# Example
```julia
analyse_output(
    "6akl",
    "/path/to/data/",
    "/path/to/data/6akl_AlphaConformer",
    "6akp"
)
```
"""
#Created the visualization plot
function analyse_output(query::String,folder_path::String,folder_model::String,objectif::String)
    #Get the dev set file

    holo_pdb = objectif*"_holo.pdb"
    apo_pdb = objectif*"_apo.pdb"
    #Found all the model predicted by AlphaFold
    model_files=find_model(folder_model::String) 

    #Compare the model to the apo and holo shape
    results=compare_model(model_files,holo_pdb,apo_pdb,folder_path,folder_model)

    #Save the CSV 
    CSV.write(folder_model*"/rmsd_results_"*query*".csv", results)
    #Visualisation of the result 
    scatter(results.RMSD_Holo,results.RMSD_Apo,
    xlabel="Holo", ylabel="Apo",
    title="Compare model with Apo and Holo form for "*query,
    marker=:circle, legend=false,
    xlims=(0,13), ylims=(0, 13))
    #Save the plot
    savefig(folder_model*"/rmsd_scatter_"*query*".png")  # Sauvegarde du plot
end

#Get input parameter from command line
function parse_commandline()
    s = ArgParseSettings(description = "Analyse AlphaConformers output: compare predicted models to apo and holo structures.")

    @add_arg_table! s begin
        "--query", "-q"
            help = "Name of the PDB entry to analyse (e.g. 6akl)"
            arg_type = String
            required = true
        "--objectif", "-o"
            help = "Basename for apo/holo reference files (e.g. 6akp)"
            arg_type = String
            required = true
        "--folder_path", "-p"
            help = "Path to the folder containing apo and holo PDB files"
            arg_type = String
            required = true
        "--folder_model", "-m"
            help = "Path to the AlphaConformers output folder (default: <folder_path>/<query>_AlphaConformer)"
            arg_type = String
            default = nothing
    end

    return parse_args(s)
end

################################################ MAIN ################################################
#Run script with command line arguments, e.g.:
#julia output_AlphaConformers.jl --query 6akl --objectif 6akp --folder_path /stockage/.../data/
args = parse_commandline()

query        = args["query"]
objectif     = args["objectif"]
folder_path  = args["folder_path"]
folder_model = something(args["folder_model"], folder_path * query * "_AlphaConformer")

analyse_output(query, folder_path, folder_model, objectif)


########################################### END ##############################################################

