# Analyse output of AlphaConformer
using MIToS.PDB
using DataFrames, CSV
using Plots
using AlphaConformers
using Glob

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
        scores_file = joinpath(cluster, "af/predictions/scores.tsv")   # Chemin vers les PDB
        println(models_path)
        if isdir(models_path) && isfile(scores_file)  # Vérifier si le chemin existe
            df_scores = CSV.read(scores_file, DataFrame, delim='\t')

            max_pLDDT = maximum(df_scores.pLDDT)
            best_models = df_scores[df_scores.pLDDT .== max_pLDDT, :]
            println(best_models)
            pdb_files = glob("sequences_unrelaxed_rank_*_alphafold2_ptm_model_*_seed_*.pdb", models_path)
            println(length(pdb_files))
            for row in eachrow(best_models)
                rank_str = lpad(row.rank, 3, '0')  # Formater le rank en "009"
                matching_pdb = filter(pdb -> occursin("rank_" * rank_str * "_", pdb), pdb_files)
                println(length(matching_pdb))
                # Ajouter les fichiers trouvés
                append!(model_files, matching_pdb)
            end
        end
    end
    return model_files
end
function analyse_output(query::String,folder_path::String)
    df_info = CSV.read(folder_path*"/info_dev_set.csv", DataFrame, delim=',')
    println(df_info)
    filtered_rows = filter(row -> occursin(query, row.PDB_apo), df_info)
    println(filtered_rows)
    row = first(filtered_rows, 1)  # Prend la première ligne

    apo_pdb = string(row.PDB_apo[1], "_", row.CHAIN_apo[1], "_", row.INDEX_apo[1], ".pdb.gz")
    holo_pdb = string(row.PDB_holo[1], "_", row.CHAIN_holo[1], "_", row.INDEX_holo[1], ".pdb.gz")

    # Construire le chemin du dossier modèle
    folder_model=folder_path*query*"/"
    println("Apo PDB Path: ", apo_pdb)
    println("Holo PDB Path: ", holo_pdb)
    println("Folder Model: ", folder_model)

    model_files=find_model(folder_model::String)
    println(length(model_files))
    println(typeof(model_files))

    results=compare_model(model_files,holo_pdb,apo_pdb,folder_path,folder_model)
    CSV.write(folder_path*"/"*query*"/rmsd_results_"*query*".csv", results)

    #Visualisation of the result 
    scatter(results.RMSD_Holo, results.RMSD_Apo,
    xlabel="Holo", ylabel="Apo",
    title="Compare model with Apo and Holo form for "*query,
    marker=:circle, legend=false,
    xlims=(0, 7), ylims=(0, 7))

    savefig(folder_path*"/"*query*"/rmsd_scatter_"*query*".png")  # Sauvegarde du plot
end

folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
query="1FMF"
analyse_output(query,folder_path)
