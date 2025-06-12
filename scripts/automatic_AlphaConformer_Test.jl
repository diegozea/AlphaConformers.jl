#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=38
#SBATCH --output=automatic_AlphaConformer_Test_results.jl.o%j.out
=#


import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")
Pkg.status("MIToS")
using MIToS
using MIToS.PDB
using DataFrames
import CSV
using Revise
using AlphaConformers
using Glob


##### Function from Analyse_output_Alpha_Conformer.jl #####
#Calculate the RMSD 
function compare_model(model_files::Vector,holo_pdb::String,apo_pdb::String,folder_path::String,folder_model::String)
    results = DataFrame(Model=String[], RMSD_Holo=Float64[], RMSD_Apo=Float64[])
    for model in model_files
        println(model)

        structure_model = read(model, PDBFile,group="ATOM")
        structure_apo = read(folder_path*apo_pdb, PDBFile,group="ATOM")
        structure_holo = read(folder_path*holo_pdb, PDBFile,group="ATOM")

        _,_,_,rmsd_holo,coverage,_=AlphaConformers.structural_alignment(structure_model, structure_holo)
        _,_,_,rmsd_apo,coverage,_ =AlphaConformers.structural_alignment(structure_model, structure_apo)
        println(basename(model))
        push!(results, (basename(model), rmsd_holo, rmsd_apo))
        
    end
    return results 
end

function find_model(folder_model::String)
    model_files = String[]
    # Trouver tous les dossiers qui commencent par "cluster"
    cluster_dirs = filter(isdir, glob("cluster*", folder_model))
    for cluster in cluster_dirs
        models_path = joinpath(cluster, "af/predictions/sequences/models")  #  Chemin des modèles
        #scores_file = joinpath(cluster, "af/predictions/scores.tsv")   # Chemin vers les PDB
        println(models_path)
        if isdir(models_path)   # Vérifier si le chemin existe
            #df_scores = CSV.read(scores_file, DataFrame, delim='\t')

            #max_pLDDT = maximum(df_scores.pLDDT)
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
function analyse_output(query::String,PATH::String,folder_model::String)
    df_info = CSV.read(PATH*"/info_dev_set.csv", DataFrame, delim=',')
    println(df_info)
    filtered_rows = filter(row -> occursin(query, row.PDB_apo), df_info)
    println(filtered_rows)
    row = first(filtered_rows, 1)  # Prend la première ligne

    apo_pdb = string(row.PDB_apo[1], "_", row.CHAIN_apo[1], "_", row.INDEX_apo[1], ".pdb.gz")
    holo_pdb = string(row.PDB_holo[1], "_", row.CHAIN_holo[1], "_", row.INDEX_holo[1], ".pdb.gz")

    # Construire le chemin du dossier modèle
    
    println("Apo PDB Path: ", apo_pdb)
    println("Holo PDB Path: ", holo_pdb)
    println("Folder Model: ", folder_model)

    model_files=find_model(folder_model::String)
    println(length(model_files))
    println(typeof(model_files))

    results=compare_model(model_files,holo_pdb,apo_pdb,PATH,folder_model)
    return results
end

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
########################################################################################################

const PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data"
cd(PATH)

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"
#const ALPHAFOLD_DB = nothing

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"

# Create folder for the test results
output = joinpath(PATH,"automatic_AlphaConformer_Test_results")

# Read the dev set file
file_path="info_dev_set.csv"
dev_set=DataFrames.DataFrame(CSV.File(file_path,
comment="#", missingstring=["", "None"])) 

# Get the dev_set information
df_info_path="Supplementary_Table_1_91_apo_holo_pairs.csv"
df_info=DataFrames.DataFrame(CSV.File(df_info_path,
        comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND

# Final CSV to save the results
csv_result= "comparaison_resultat_AlphaConformer.csv"

#The hyperparameters for the test
database = [[FOLDSEEK_DB],[FOLDSEEK_DB,ALPHAFOLD_DB]]
evalue=[1e-5, 1e-1, nothing]
cutoff=[0.5,1,1.5,2]

#For each pdb in the dev_set
for row in eachrow(dev_set)
    query = row.PDB_apo
    chain= row.CHAIN_apo
    apo_model=row.INDEX_apo
    REF_PDB = joinpath(PATH, query*"_"*chain*"_"*string(apo_model)*".pdb.gz")
    #Change the hyperparameters
    for database_path in database
        evalue_count=0
        for eval in evalue
            evalue_count += 1
            if isdir(output)
                rm(output, recursive=true)
            end
            mkdir(output)
            cutoff_count=0
            for cut in cutoff
                cutoff_count += 1
                @info "Processing query: $query, database : $database_path, evalue: $eval, cutoff: $cut"
                # Create a folder for the current query
                query_folder = joinpath(output, "$(query)_$(chain)_$(apo_model)_database_$(length(database_path))_eval_$(eval)_cutoff_$(cut)")
                if isdir(query_folder)
                    rm(query_folder; recursive=true, force=true)
                end
                mkdir(query_folder)
                # Run AlphaConformers
                try
                    @show database_path
                    AlphaConformers.alphaconformers(REF_PDB, PDB_FOLDER, query_folder; db=database_path, evalue_cutoff=eval, cutoff=cut)
                    AlphaConformers.run_alphafold(query_folder, colabfold_path=COLABFOLD_PATH)
                catch e
                    @error "Error processing $query with chain $chain and apo model $apo_model: $e"
                    continue
                end
            end
            folders_result = filter(isdir, readdir(output, join=true))
            # Attendre que tous les jobs alphafold soient terminés
            while true
                # 1. Prendre le dernier dossier de folders_result (par ordre alphabétique)
                last_folder = sort(folders_result)[end]
                @show last_folder
                # 2. Chercher le dernier dossier cluster_* dans ce dossier
                cluster_dirs = filter(isdir, sort(glob("cluster_*", last_folder)))
                last_cluster = cluster_dirs[end]
                @show last_cluster
                # 3. Aller dans af/predictions/
                pred_dir = joinpath(last_cluster, "af", "predictions")
                pred_contents = readdir(pred_dir)
                if !isempty(pred_contents) #si le dossier n'est pas vide
                    # 4. Chercher le dossier models
                    models_dir = joinpath(pred_dir, "sequences/models")
                    model_files = filter(f -> isfile(joinpath(models_dir, f)), readdir(models_dir))
                    if length(model_files) == 10
                        @info "10 modèles trouvés, on sort du while"
                        break
                    end
                end
                @info "En attente de la fin des jobs alphafold pour la requête $query, database : $database_path, evalue: $eval, cutoff: $cut"
                sleep(600) # attend 10 min avant de re-vérifier
            end
            @info "Tous les jobs alphafold sont terminés."
            # Analyse des résultats
            @info "Analyse des résultats pour la requête $query"
            
            #get the csv file
            if isfile(csv_result)
                df_csv_result=DataFrames.DataFrame(CSV.File(csv_result,
                    comment="#", missingstring=["", "None"])) 
            else
                # Create a DataFrame to store the results
                df_csv_result = DataFrames.DataFrame(Query=String[],Method = String[], Database = Int64[], Evalue = Float64[], Cutoff = Float64[], Nombre_cluster = Float64[], Mean_MSA= Float64[], Sequence_in_msa = Int128[], rmsd_apo_holo = Float64[], min_apo=Float64[], mean_apo =Float64[], var_apo=Float64[],min_holo=Float64[], mean_holo = Float64[], var_holo =Float64[])
            end

            row_match = filter(row -> occursin(query, row.apo_id), df_info)
            rmsd_apo_holo=row_match.rmsd_apo_holo[1]
            #Save the information in the csv file
            for folder in folders_result
                result_csv_AF=analyse_output(String(query),PATH,folder)
                result_directory=analyze_directory(folder,false)
                result_directory_seqs = Set(lowercase.(result_directory.unique_sequence_names))
                push!(df_csv_result,(query,"Foldseek", length(database_path),eval,cut,result_directory.n_clusters,result_directory.avg_lines_per_a3m,length(result_directory_seqs), rmsd_apo_holo,minimum(result_csv_AF.RMSD_Apo),mean(result_csv_AF.RMSD_Apo),var(result_csv_AF.RMSD_Apo),minimum(result_csv_AF.RMSD_Holo),mean(result_csv_AF.RMSD_Holo),var(result_csv_AF.RMSD_Holo)))
            end
            CSV.write(csv_result,df_csv_result) 
            break 
        end
        break
    end
    break

end
@info "Test terminé. Résultats enregistrés dans $csv_result"