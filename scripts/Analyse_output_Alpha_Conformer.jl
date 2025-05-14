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
function analyse_output(query::String,folder_path::String)
    df_info = CSV.read(folder_path*"/info_dev_set.csv", DataFrame, delim=',')
    println(df_info)
    filtered_rows = filter(row -> occursin(query, row.PDB_apo), df_info)
    println(filtered_rows)
    row = first(filtered_rows, 1)  # Prend la première ligne

    apo_pdb = string(row.PDB_apo[1], "_", row.CHAIN_apo[1], "_", row.INDEX_apo[1], ".pdb.gz")
    holo_pdb = string(row.PDB_holo[1], "_", row.CHAIN_holo[1], "_", row.INDEX_holo[1], ".pdb.gz")

    # Construire le chemin du dossier modèle
    folder_model=folder_path*query*"_Update"
    println("Apo PDB Path: ", apo_pdb)
    println("Holo PDB Path: ", holo_pdb)
    println("Folder Model: ", folder_model)

    model_files=find_model(folder_model::String)
    println(length(model_files))
    println(typeof(model_files))

    results=compare_model(model_files,holo_pdb,apo_pdb,folder_path,folder_model)
    CSV.write(folder_path*"/"*query*"_Update/rmsd_results_"*query*".csv", results)

    #Visualisation of the result 
    scatter(results.RMSD_Apo, results.RMSD_Holo,
    xlabel="Apo", ylabel="Holo",
    title="Compare model with Apo and Holo form for "*query,
    marker=:circle, legend=false,
    xlims=(0,7), ylims=(0, 6))

    savefig(folder_path*"/"*query*"_Update/rmsd_scatter_"*query*".png")  # Sauvegarde du plot
end
using Glob
using Statistics

function analyze_directory(base_dir::String,div::Bool)
    cluster_dirs = filter(isdir, sort(glob("cluster_*", base_dir)))
    n_clusters = length(cluster_dirs)
    template_counts = Int[]
    a3m_line_counts = Int[]
    template_files_set = Set{String}()

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
            lines = countlines(a3m_file)
            push!(a3m_line_counts, lines)
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
        template_files = template_files_set
    )
end

function compare_directories(dir1::String, dir2::String)
    stats1 = analyze_directory(dir1,true)
    stats2 = analyze_directory(dir2,false)

    println("=== Comparaison de $dir1 et $dir2 ===")
    println("Clusters: $(stats1.n_clusters) vs $(stats2.n_clusters)")
    println("Fichiers templates moyens: $(stats1.avg_templates_per_cluster) vs $(stats2.avg_templates_per_cluster)")
    println("Lignes moyennes .a3m: $(stats1.avg_lines_per_a3m) vs $(stats2.avg_lines_per_a3m)")

    files1_lower = Set(lowercase.(collect(stats1.template_files)))
    files2_lower = Set(lowercase.(collect(stats2.template_files)))

    only_in_dir1 = setdiff(files1_lower, files2_lower)
    only_in_dir2 = setdiff(files2_lower, files1_lower)

    println("\nFichiers uniquement dans $dir1/templates:")
    println(length(only_in_dir1))
    println(length(files1_lower))
    println("\nFichiers uniquement dans $dir2/templates:")
    println(length(only_in_dir2))
    println(length(files2_lower))
end

folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
query="1AKZ"
#analyse_output(query,folder_path)

chemin_dossier_1 = joinpath(folder_path,"1AKZ")
chemin_dossier_2 = joinpath(folder_path,"1AKZ_Update")
resultats = compare_directories(chemin_dossier_1, chemin_dossier_2)