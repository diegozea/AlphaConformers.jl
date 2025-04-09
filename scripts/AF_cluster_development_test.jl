#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#PBS -l host=node48
#PBS -l walltime=900:00:00
#PBS -l mem=100gb
#PBS -l ncpus=40
#PBS -j oe
=#
import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
Pkg.add("Glob")
Pkg.add("DataFrames")
using Glob
using DataFrames, CSV
using Plots
using MIToS.PDB
using AlphaConformers

function AF_cluster(query,create_file)
    ### Create the folder to save the result
    folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
    cluster_dir=folder_path*query*"_AF_CLUSTER_template"
    if create_file
        if !isdir(cluster_dir)
            mkdir(cluster_dir)
        end

        msa=joinpath(cluster_dir, "msa.a3m")
        if !isfile(msa)
            msa_origin=folder_path*"/"*query*"/afdb_up_results/msa.a3m"
            cp(msa_origin, msa)
        end

        PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/"
        cd(PATH)

        ### Create the cluster 
        run(`python3 scripts/ClusterMSA.py EX -i $msa -o msas`)
        msas=joinpath(cluster_dir, "msas")
        if !isdir(msas)
            cp("msas", msas)
        end

        ### Run AF2
        #joinpath
        output_dir=joinpath(cluster_dir,"output")
        template=cluster_dir*"/template/"
        isdir(template) || mkdir(template)

        ## Move the template from AlphaConformer
        cluster_dirs = filter(isdir, glob("cluster*", folder_path*query))
        for cluster in cluster_dirs
            files = glob("*", cluster*"/templates")  # Tous les fichiers du sous-dossier
            for file in files
                cp(file, joinpath(template, basename(file)), force=true)
            end
        end
        println("✅ Copie terminée !")
        files = glob("*", template)

        for file in files
            if !occursin(r"^[a-zA-Z0-9]{4}\.(pdb|cif)$", basename(file))
                rm(file)  # Supprime le fichier
            end
        end
        println("✅ Nettoyage terminé !")
    else 
        if !isdir(cluster_dir)
            println("ERROR : the folder $cluster_dir doesn't exist ")
            return nothing
        end
        msas=joinpath(cluster_dir, "msas")
        output_dir=joinpath(cluster_dir,"output")
        template=cluster_dir*"/template/"
    end

    COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
    #af_command = `$COLABFOLD_PATH $msas $output_dir --msa-input`
    af_command = `$COLABFOLD_PATH $msas $output_dir --use-templates 1 --msa-input --custom-template-path $template --overwrite-existing-results`
    @info "Running AlphaFold command: $af_command"
    run(af_command)
end

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

function visualisation(folder_model,holo_pdb,apo_pdb,folder_path,query,template)
    model_files = String[]
    # Trouver tous les dossiers qui commencent par "cluster"
    cluster_dirs = filter(f -> occursin(r"^EX_(?!U100)", basename(f)), glob("EX_U10*",folder_model))
    for cluster in cluster_dirs
        models_path = joinpath(cluster, "models")  # 📂 Chemin des modèles
        println(models_path)
        if isdir(models_path)   # Vérifier si le chemin existe
            matching_pdb = filter(isfile, glob("*.pdb", models_path))
            append!(model_files, matching_pdb)
        end
    end

    println(length(model_files))
    println(typeof(model_files))
    results=compare_model(model_files,holo_pdb,apo_pdb,folder_path,folder_model)
    CSV.write(folder_path*"/"*query*"_AF_CLUSTER/rmsd_results_"*query*"_U10.csv", results)

    #Visualisation of the result 
    scatter(results.RMSD_Holo, results.RMSD_Apo,
        xlabel="Holo", ylabel="Apo",
        title="Compare model with Apo and Holo form for $query with AF-CLUSTER U10",
        marker=:circle, legend=false,
        xlims=(0, 7), ylims=(0, 7))

    savefig(folder_path*"/"*query*"_AF_CLUSTER/rmsd_scatter_"*query*"_U10.png")  # Sauvegarde du plot

    model_files = String[]
    # Trouver tous les dossiers qui commencent par "cluster"
    cluster_dirs = filter(isdir, glob("EX_U100*", folder_model))
    for cluster in cluster_dirs
        models_path = joinpath(cluster, "models")  # 📂 Chemin des modèles
        println(models_path)
        if isdir(models_path)   # Vérifier si le chemin existe
            matching_pdb = filter(isfile, glob("*.pdb", models_path))
            append!(model_files, matching_pdb)
        end
    end

    println(length(model_files))
    println(typeof(model_files))
    results=compare_model(model_files,holo_pdb,apo_pdb,folder_path,folder_model)
    CSV.write(folder_path*"/"*query*"_AF_CLUSTER/rmsd_results_"*query*"_U100.csv", results)

    #Visualisation of the result 
    scatter(results.RMSD_Holo, results.RMSD_Apo,
        xlabel="Holo", ylabel="Apo",
        title="Compare model with Apo and Holo form for $query with AF-CLUSTER U100",
        marker=:circle, legend=false,
        xlims=(0, 7), ylims=(0, 7))

    savefig(folder_path*"/"*query*"_AF_CLUSTER/rmsd_scatter_"*query*"_U100.png")  # Sauvegarde du plot

    model_files = String[]
    cluster_labels = String[]  # Stocke le cluster d'origine pour chaque fichier
    cluster_dirs = filter(f -> occursin(r"^EX_(?!U)", basename(f)), glob("EX_*", folder_model))

    for cluster in cluster_dirs
        models_path = joinpath(cluster, "models")  # 📂 Chemin des modèles
        println(models_path)
        if isdir(models_path)   # Vérifier si le chemin existe
            matching_pdb = filter(isfile, glob("*.pdb", models_path))
            append!(model_files, matching_pdb)
            append!(cluster_labels, fill(basename(cluster), length(matching_pdb)))  # Associe chaque fichier à son cluster
        end
    end

    println(length(model_files))
    println(typeof(model_files))

    # Comparaison des modèles et récupération des résultats
    results = compare_model(model_files, holo_pdb, apo_pdb, folder_path, folder_model)
    CSV.write(folder_path * "/"*query*"_AF_CLUSTER/rmsd_results_"*query*".csv", results)

    # Attribution d'une couleur unique par cluster
    unique_clusters = unique(cluster_labels)
    color_map = Dict(cluster => getindex(distinguishable_colors(length(unique_clusters)), i) for (i, cluster) in enumerate(unique_clusters))
    colors = [color_map[cluster] for cluster in cluster_labels]  # Associe la couleur à chaque point

    # Visualisation du résultat avec couleurs spécifiques
    scatter(results.RMSD_Holo, results.RMSD_Apo,
        xlabel="Holo", ylabel="Apo",
        title="Compare model with Apo and Holo form for $query with AF-CLUSTER",
        marker=:circle, legend=true,
        xlims=(0, 7), ylims=(0, 7),
        color=colors, label=cluster_labels)

    savefig(folder_path * "/"*query*"_AF_CLUSTER/rmsd_scatter_"*query*".png")  # Sauvegarde du plot
end
################### MAIN #######################

folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
query="1AKZ"
create_file=false
## Alpha-Cluster
AF_cluster(query,create_file)
"""
#Visulaize the result 
df_info = CSV.read(folder_path*"/info_dev_set.csv", DataFrame, delim=',')
println(df_info)
filtered_rows = filter(row -> occursin(query, row.PDB_apo), df_info)
println(filtered_rows)
row = first(filtered_rows, 1)  # Prend la première ligne

apo_pdb = string(row.PDB_apo[1], "_", row.CHAIN_apo[1], "_", row.INDEX_apo[1], ".pdb.gz")
holo_pdb = string(row.PDB_holo[1], "_", row.CHAIN_holo[1], "_", row.INDEX_holo[1], ".pdb.gz")
#1AKZ,A,1,1SSP,E,1
template=true

folder_model=folder_path*query*"_AF_CLUSTER_template/output/predictions/"

visualisation(folder_model,holo_pdb,apo_pdb,folder_path,query,template)
"""