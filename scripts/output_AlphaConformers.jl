#!/stockage/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#SBATCH --nodelist=node48
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2
#SBATCH --output=output_AlphaConformers.jl.o%j.out
=#
import Pkg
Pkg.activate("/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")


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
Compare the model predicted to the apo and holo shape 
Use find_model and compare_model function 
Create a CSV and a scatter plot that is save in the output folder of ALphaConformer
"""
#Created the visualization plot
function analyse_output(query::String,folder_path::String,folder_model::String,df_info)
    #Get the dev set file
    
    filtered_rows = filter(row -> occursin(query, row.PDB_apo), df_info)
    row = first(filtered_rows, 1)  # Take the first row

    #Extract the apo and holo shape
    apo_pdb = string(row.PDB_apo[1], "_", row.CHAIN_apo[1],"_", row.INDEX_apo[1], ".pdb.gz")
    holo_pdb = string(row.PDB_holo[1], "_", row.CHAIN_holo[1],"_", row.INDEX_holo[1], ".pdb.gz")

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


################################################ MAIN ################################################
"""
Analyse the output of AlphaConformers 
Can be use to create he vsualization ad the csv to compre model predicted to apo and holo shape 
Or can be use to compare the clustering created by AlphaConformer for two different folder

Input : 
- folder_path: Path to the folder where apo and holo file are save 
- folder_model : Folder output path from ALphaConformers
- file ou query : name of the pdb to analyse 
Output : 
Plot and CSV save in the ALphaConformers folder 
and/or
Number of cluster, ligne in MSA, template, and clustering for 2 AlphaConformer result 

Analyse all our result with analyse_output function 
Need to change the plot to color by pLLDT 
"""
########################## Information to fill #################################
#Path to the folder where apo and holo file are save 
folder_path ="/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
folder_model = "/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"*query*"_AlphaConformer"
#OR file name 
query="6akl"
################################################################################

#For one run 
#Get the comparaison between every model and apo and hlo shape 
analyse_output(query,folder_path,folder_model)

########################################### END ##############################################################

