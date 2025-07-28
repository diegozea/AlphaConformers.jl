#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --output=automatic_AlphaConformer_Test_results.jl.o%j.out
=#


import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")
Pkg.status("MIToS")

# Load necessary packages 
using MIToS
using MIToS.PDB
using DataFrames
import CSV
using Revise
using AlphaConformers
using Glob
using Statistics


############################### Function from Analyse_output_Alpha_Conformer.jl #######################################

"""
Compare all the model predicted by AlphaFold to the apo and holo shape
Use the function structural_alignment from AlphaConformers
"""
#Calculate the RMSD 
function compare_model(model_files::Vector,holo_pdb::String,apo_pdb::String,folder_path::String,folder_model::String)
    results = DataFrame(Model=String[], RMSD_Holo=Float64[], RMSD_Apo=Float64[])
    for model in model_files

        structure_model = read_file(model, PDBFile,group="ATOM")
        structure_apo = read_file(folder_path*apo_pdb, PDBFile,group="ATOM")
        structure_holo = read_file(folder_path*holo_pdb, PDBFile,group="ATOM")

        _,_,_,rmsd_holo,coverage,_=AlphaConformers.structural_alignment(structure_model, structure_holo)
        _,_,_,rmsd_apo,coverage,_ =AlphaConformers.structural_alignment(structure_model, structure_apo)
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
    plddt=Float64[]
    # Trouver tous les dossiers qui commencent par "cluster"
    cluster_dirs = filter(isdir, glob("cluster*", folder_model))
    for cluster in cluster_dirs
        models_path = joinpath(cluster, "af/predictions/sequences/models")  #  Chemin des modèles
        scores_file = joinpath(cluster, "af/predictions/scores.tsv")   # Chemin vers les PDB
        if isdir(models_path)   # Vérifier si le chemin existe
            df_scores = CSV.read(scores_file, DataFrame, delim='\t')

            max_pLDDT = maximum(df_scores.pLDDT)
            push!(plddt,max_pLDDT)
            #best_models = df_scores[df_scores.pLDDT .== max_pLDDT, :]
            #println(best_models)
            pdb_files = glob("sequences_unrelaxed_rank_*_alphafold2_ptm_model_*_seed_*.pdb", models_path)
            append!(model_files, pdb_files)
        else 
            println("Path not found")
        end
    end
    return model_files, plddt
end

"""
Compare the model predicted to the apo and holo shape 
Use find_model and compare_model function 
Create a CSV and a scatter plot that is save in the output folder of ALphaConformer
"""
#Created the visualization plot
function analyse_output(query::String,PATH::String,folder_model::String)
    df_info = CSV.read(PATH*"/info_dev_set.csv", DataFrame, delim=',')
    filtered_rows = filter(row -> occursin(query, row.PDB_apo), df_info)
    
    row = first(filtered_rows, 1)  # Prend la première ligne

    apo_pdb = string(row.PDB_apo[1], "_", row.CHAIN_apo[1], "_", row.INDEX_apo[1], ".pdb.gz")
    holo_pdb = string(row.PDB_holo[1], "_", row.CHAIN_holo[1], "_", row.INDEX_holo[1], ".pdb.gz")

    # Construire le chemin du dossier modèle
    

    model_files,plddt=find_model(folder_model::String)
    

    results=compare_model(model_files,holo_pdb,apo_pdb,PATH,folder_model)
    return results,plddt
end

"""
Analyse the output of AlphaConformer
Register the number of cluster, the number of file use , the lenght of the MSA...
"""
#Analyse the output directory of AlphaConformer 
function analyze_directory(base_dir::String,div::Bool)
    cluster_dirs = filter(isdir, sort(glob("cluster_*", base_dir)))
    n_clusters = length(cluster_dirs)
    template_counts = Int[]
    a3m_line_counts = Int[]
    template_files_set = Set{String}()
    all_seq_names = Set{String}() 

    for cluster_dir in cluster_dirs
        template_dir = joinpath(cluster_dir, "templates")
        a3m_files = glob("*.a3m", cluster_dir)
        
        # Nombre de fichiers dans templates
        if isdir(template_dir)
            templates = readdir(template_dir)
            push!(template_counts, length(templates))
            foreach(f -> push!(template_files_set, f), templates)
        else
            push!(template_counts, 0)
        end

        # Nombre de lignes dans les fichiers .a3m
        for a3m_file in a3m_files
            msa=MIToS.MSA.read_file(a3m_file, MIToS.MSA.FASTA, generatemapping=true)
            seqnames = MIToS.MSA.sequencenames(msa)[2:end]
            push!(a3m_line_counts, length(seqnames))
            union!(all_seq_names, seqnames)
        end
    end
    if div 
        avg_templates = mean(template_counts)/2
    else
        avg_templates = mean(template_counts)
    end
    avg_a3m_lines = isempty(a3m_line_counts) ? 0 : mean(a3m_line_counts)

    return (
        n_clusters = n_clusters,
        avg_templates_per_cluster = avg_templates,
        avg_lines_per_a3m = avg_a3m_lines,
        template_files = template_files_set,
        unique_sequence_names = collect(all_seq_names)
    )
end
######################################### MAIN ##############################################
"""
This code execute AlphaConformers with AF2 for different parameter 
Three parameter can be test : database, evalue and cutoff
Do the analysis of the output : number of cluster, number of sequence, lenght MSA

Input : 
- PATH : Path to the folder where apo and holo file are save 
- PDB_FOLDER : Define the path to the PDB folder
- FOLDSEEK_DB: Define the path to the foldseek database folder
- ALPHAFOLD_DB : Define the path to the afdb folder
- COLABFOLD_PATH : Path to the ColabFold script
- output : Create folder for the results
- database, evalue, cutoff : The hyperparameters for the test
- csv_result : Final CSV to save the results
Output :
CSV that registrer for each parameter the performance of AlphaConformers

Create the folder for each cutoff than delete the folder for the next evalue and/or database test
Use less stockage space 

Very long to execute 
Depend directly of the number of cluster 
"""
########################## Information to fill #################################
#Path to the folder where apo and holo file are save 
const PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
cd(PATH)
# Define the path to the PDB folder
const PDB_FOLDER = "/alpha/database/pdb/pdb_files"
# Define the path to the foldseek database folder
const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"
# Define the path to the afdb folder
const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"
# Path to the ColabFold script
const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
# Create folder for the results
output = joinpath(PATH,"automatic_AlphaConformer_Test_results")
                ##################################
# Read the dev set file
file_path="info_dev_set.csv" #Created in  AlphaConformers_development_test.jl
dev_set=DataFrames.DataFrame(CSV.File(file_path,
comment="#", missingstring=["", "None"])) 

# Get the dev_set information
df_info_path="Supplementary_Table_1_91_apo_holo_pairs.csv"
df_info=DataFrames.DataFrame(CSV.File(df_info_path,
        comment="#", missingstring=["", "None"])) 

#The hyperparameters for the test
database = [[FOLDSEEK_DB],[FOLDSEEK_DB,ALPHAFOLD_DB]]
evalue=[1e-5, NaN]
cutoff=[1,1.5,2]

# Final CSV to save the results
csv_result= "comparaison_resultat_AlphaConformer.csv"
##################################################################################
@show "Start"

#get the csv file
if isfile(csv_result)
    df_csv_result=DataFrames.DataFrame(CSV.File(csv_result,
        comment="#", missingstring=["", "None"])) 
else
    # Create a DataFrame to store the results
    df_csv_result = DataFrames.DataFrame(Query=String[],Method = String[], Database = Int64[], Evalue = Union{Float64, Missing}[], Cutoff = Float64[], Nombre_cluster = Float64[], Mean_MSA= Float64[], Sequence_in_msa = Int128[], rmsd_apo_holo = Float64[], mean_plddt = Union{Missing,Float64}[], min_apo=Float64[], mean_apo =Float64[], var_apo=Float64[],min_holo=Float64[], mean_holo = Float64[], var_holo =Float64[])
end


#For each pdb in the dev_set
for row in eachrow(dev_set)
    query = row.PDB_apo
    if query =="6OY9" || query=="2UZ5"
        continue # Skip those query
    end
    chain= row.CHAIN_apo
    apo_model=row.INDEX_apo
    REF_PDB = joinpath(PATH, query*"_"*chain*"_"*string(apo_model)*".pdb.gz")
    #Change the hyperparameters
    for database_path in database
        
        for eval in evalue
            # Check the test already done 
            already_done = any(row -> row[1] == query &&
                                    row[2] == "Foldseek" &&
                                    row[3] == length(database_path) &&
                                    ((isnan(row[4]) && isnan(eval)) || row[4] == eval), eachrow(df_csv_result))
            if already_done
                @info "Skipping query: $query, database : $database_path, evalue: $eval - already processed."
                continue
            end

            if isdir(output)
                rm(output, recursive=true) #keep only the folder for the new test
            end
            mkdir(output)
            
            for cut in cutoff
                #Same code as in AlphaConformers_development_test.jl
                @info "Processing query: $query, database : $database_path, evalue: $eval, cutoff: $cut"
                # Create a folder for the current query
                query_folder = joinpath(output, "$(query)_$(chain)_$(apo_model)_database_$(length(database_path))_eval_$(eval)_cutoff_$(cut)")
                if isdir(query_folder)
                    continue # Skip if the folder already exists
                    rm(query_folder; recursive=true, force=true)
                end
                mkdir(query_folder)
                # Run AlphaConformers
                try 
                    AlphaConformers.alphaconformers(REF_PDB, PDB_FOLDER, query_folder; db=database_path, evalue_cutoff=eval, cutoff=cut)
                    AlphaConformers.run_alphafold(query_folder, colabfold_path=COLABFOLD_PATH)
                catch e
                    @error "Error processing query: $query, database : $database_path, evalue: $eval, cutoff: $cut - $(e)"
                    continue # Skip to the next query   
                end
                
            end
            
            folders_result = filter(isdir, readdir(output, join=true))
            # Wait until AF2 is finish
            while true
                # 1. Check the last folder
                last_folder = sort(folders_result)[end]
                # 2. Check the last cluster file
                cluster_dirs = filter(isdir, sort(glob("cluster_*", last_folder)))
                sorted_clusters = sort(cluster_dirs, by = x -> parse(Int, split(x, "_")[end]))
                last_cluster = sorted_clusters[end]
                # 3. Go in  af/predictions/
                pred_dir = joinpath(last_cluster, "af", "predictions")
                pred_contents = readdir(pred_dir)
                if !isempty(pred_contents) #if the folder exist 
                    # 4. Get the model
                    models_dir=joinpath(pred_dir, "sequences","models")
                    if isdir(models_dir)
                        model_files = filter(f -> isfile(joinpath(models_dir, f)), readdir(models_dir))
                        if length(model_files) == 10 #check that we hav the 10 predicted
                            break
                        end
                    else
                        @warn "models_dir does not exist: $models_dir"
                        sleep(600) #Job not finished we wait
                        continue
                    end
                end
                @info "Wait until the end of AlphaFold for $query, database : $database_path, evalue: $eval"
                sleep(600) # Job not finished we wait
            end
            @info "AF2 have finish"

            # Result analysis - 

            row_match = filter(row -> occursin(query, row.apo_id), df_info)
            rmsd_apo_holo=row_match.rmsd_apo_holo[1]
            #Save the information in the csv file
            index=0
            #check every output folder of ALphaConformer
            for folder in folders_result
                index += 1
                result_csv_AF,plddt=analyse_output(String(query),PATH,folder)
                
                if isempty(plddt)
                    @warn "No pLDDT values found for query $query in folder $folder"
                    push!(df_csv_result,(query,"Foldseek", length(database_path),eval,cutoff[index],result_directory.n_clusters,result_directory.avg_lines_per_a3m,length(result_directory_seqs), rmsd_apo_holo,missing,missing,missing,missing,missing,missing,missing))
                    continue
                end
                #code use in Analyse_output_Alpha_Conformer.jl
                result_directory=analyze_directory(folder,false)
                result_directory_seqs = Set(lowercase.(result_directory.unique_sequence_names))
                push!(df_csv_result,(query,"Foldseek", length(database_path),eval,cutoff[index],result_directory.n_clusters,result_directory.avg_lines_per_a3m,length(result_directory_seqs), rmsd_apo_holo,mean(plddt),minimum(result_csv_AF.RMSD_Apo),mean(result_csv_AF.RMSD_Apo),std(result_csv_AF.RMSD_Apo),minimum(result_csv_AF.RMSD_Holo),mean(result_csv_AF.RMSD_Holo),std(result_csv_AF.RMSD_Holo)))
            end
            CSV.write(csv_result,df_csv_result) #save the comparaison values in  the csv 
             # On ne traite qu'une requête à la fois
        end
        
    end
    
end
@info "Test finish, Result register in  $csv_result"
@show "End"
############################################################################################################################