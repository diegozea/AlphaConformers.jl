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

"""
    run_alphafold(clusters_folder::String; colabfold_path::String=get(ENV, "COLABFOLD_PATH", ""))

This function runs AlphaFold 2 (colabfold_batch) for each cluster in the 
`clusters_folder`. The `colabfold_path` argument is the path to the colabfold_batch program.
Alternatively, you can set the COLABFOLD_PATH environment variable.
"""
#Function that Run colabfold for every folder 
function run_alphafold(clusters_folder::String; colabfold_path::String=get(ENV, "COLABFOLD_PATH", ""))
    if isempty(colabfold_path)
        throw(ErrorException("The path to ColabFold is not defined, please set the COLABFOLD_PATH environment variable or the colabfold_path keyword argument."))
    end
    cluster_folders = filter!(
        dir -> occursin("EX_", dir), #For the output of AF-Cluster
        readdir(clusters_folder, join=true))
    if isempty(cluster_folders)
        throw(ErrorException("No EX_* folders were found in $clusters_folder"))
    end
    # remember the current working directory
    current_dir = pwd()
    try
        # run AlphaFold for each cluster
        for folder in cluster_folders
            cd(folder)
            msa=basename(folder)
            println(msa)
            @info "Running AlphaFold for $(abspath(folder))"
            if isdir("template")
                if !isempty(readdir("template"))
                    af_command = `$colabfold_path $msa af --use-templates 1 --msa-input --custom-template-path template/ --num-seeds 5 --use-dropout --num-models 2 --overwrite-existing-results`
                    @info "Running AlphaFold command: $af_command"
                    run(af_command)
                else
                    @warn "No templates found in $(abspath("template"))"
                end
            else
                @warn "There is no templates folder in $(abspath(folder))"
            end
        end
    finally
        # return to the original working directory
        cd(current_dir)
    end
end
"""
Run AF-Cluster to create the MSA and run the colabfold
Take in input the name of the query and a boolean to know if the folder is alrready there 
"""
#Create the folder and Run colabfold
function AF_cluster(query,create_file)
    ### Create the folder to save the result
    folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
    folder_dir=folder_path*query*"_AF_CLUSTER"
    if create_file
        if !isdir(folder_dir)
            mkdir(folder_dir)
        else 
            rm(expanduser(folder_dir); force=true, recursive=true)
            mkdir(folder_dir)
        end

        #take the msa use in AlphaConformers
        msa=joinpath(folder_path, "2QKEE_colabfold.a3m")
        if !isfile(msa)
            msa_origin=folder_path*"/"*query*"/afdb_up_results/msa.a3m"
            cp(msa_origin, msa)
        end

        PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/"
        cd(PATH)

        ### Create the cluster with AF-cluster
        run(`python3 scripts/ClusterMSA.py EX -i $msa -o msas`)
        msas=joinpath(folder_dir, "msas")
        if !isdir(msas)
            cp("msas", msas)
        end

        ### Order the folder 
        #Create the output folder
        output_dir=folder_dir*"/output"
    else 
        #if we don't create the file make sure the folder exist
        if isdir(folder_dir)
            println("ERROR : the folder $folder_dir doesn't exist ")
            return nothing
        end
        output_dir=joinpath(folder_dir,"output")
        msas=joinpath(folder_dir, "msas")

    end
    ##Run AF2
    COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
    af_command = `$COLABFOLD_PATH $msas $output_dir --msa-input`
    @info "Running AlphaFold command: $af_command"
    run(af_command)
    
end
"""
Run AF-Cluster to create the MSA and order the file to use run_alphafold()
Take in input the name of the query and a boolean to know if the folder is alrready there 
"""
#Create the folder and Run colabfold
function AF_cluster_template(query,create_file)
    ### Create the folder to save the result
    folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
    folder_dir=folder_path*query*"_AF_CLUSTER_template"
    if create_file
        if !isdir(folder_dir)
            mkdir(folder_dir)
        else 
            rm(expanduser(folder_dir); force=true, recursive=true)
            mkdir(folder_dir)
        end

        #take the msa use in AlphaConformers
        msa=joinpath(folder_dir, "msa.a3m")
        if !isfile(msa)
            msa_origin=folder_path*"/"*query*"/afdb_up_results/msa.a3m"
            cp(msa_origin, msa)
        end

        PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/"
        cd(PATH)

        ### Create the cluster with AF-cluster
        run(`python3 scripts/ClusterMSA.py EX -i $msa -o msas`)
        msas=joinpath(folder_dir, "msas")
        if !isdir(msas)
            cp("msas", msas)
        end

        ### Order the folder 
        #Create the output folder
        output_dir=joinpath(folder_dir,"output")
        isdir(output_dir) || mkdir(output_dir)
        cluster_dirs = filter(isfile, glob("*a3m", msas))

        ## Get the template from AlphaConformer
        template_file=[]
        cluster_dirs_AC = filter(isdir, glob("cluster*", folder_path*query))
        for cluster in cluster_dirs_AC
            files = glob("*", cluster*"/templates")  # Tous les fichiers du sous-dossier
            for file in files
                if occursin(r"^[a-zA-Z0-9]{4}\.(pdb|cif)$", basename(file))
                    push!(template_file,file)
                end
            end
        end
        println(length(template_file))
        #For every msa created by AF-Cluster 
        for cluster in cluster_dirs
            ##Create the folder as in AlphaConformers
            #Create the folder for each msa 
            output_dir_step=joinpath(output_dir,basename(cluster))
            isdir(output_dir_step) || mkdir(output_dir_step)
            println(output_dir_step)
            #paste the .a3m file in the folder 
            cp(joinpath(msas, basename(cluster)),joinpath(output_dir_step, basename(cluster)), force=true)
            #Create the folder af to take the output of colabfold
            output_dir_step_af=joinpath(output_dir_step,"af")
            isdir(output_dir_step_af) || mkdir(output_dir_step_af)
            #Create the folder to save the template
            template=joinpath(output_dir_step,"template")
            isdir(template) || mkdir(template)
            #paste the template find in AlphaConformer 
            for file in template_file
                cp(file, joinpath(template, basename(file)), force=true) #in upper case 
                cp(file, joinpath(template, lowercase(basename(file))), force=true) #in lower case 
            end
            println("✅ Copie terminée !")
        end
    else 
        #if we don't create the file make sure the folder exist
        if !isdir(folder_dir)
            println("ERROR : the folder $folder_dir doesn't exist ")
            return nothing
        end
        output_dir=joinpath(folder_dir,"output")
    end
    ## Run AF2
    COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
    run_alphafold(output_dir,colabfold_path=COLABFOLD_PATH)
    
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
        if template 
            models_path = cluster * "/af/predictions/" * split(basename(cluster), ".")[1] * "/models"  # 📂 Chemin des modèles
        else 
            models_path = joinpath(cluster, "models")
        end
        println(models_path)
        if isdir(models_path)   # Vérifier si le chemin existe
            matching_pdb = filter(isfile, glob("*.pdb", models_path))
            append!(model_files, matching_pdb)
        end
    end

    println(length(model_files))
    println(typeof(model_files))
    results=compare_model(model_files,holo_pdb,apo_pdb,folder_path,folder_model)
    if template 
        CSV.write(folder_path*"/"*query*"_AF_CLUSTER_template/rmsd_results_"*query*"_U10.csv", results)
    else 
        CSV.write(folder_path*"/"*query*"_AF_CLUSTER/rmsd_results_"*query*"_U10.csv", results)
    end

   scatter(
        results.RMSD_Apo,
        results.RMSD_Holo,
        xlabel = "Apo",
        ylabel = "Holo",
        title = "$query with AF-CLUSTER U10",
        marker = :circle,
        legend = false,
        xlims = (0, 30),
        ylims = (0, 30)
    )

    # Création de l'inset (zoom sur 0–7)
    scatter!(
        results.RMSD_Apo,
        results.RMSD_Holo,
        xlims = (0, 7),
        ylims = (0, 7),
        inset = (1, bbox(0.05, 0.05, 0.5, 0.5, :bottom, :right)),
        subplot = 2,
        marker = :circle,
        markersize = 3,  
        legend = false,
        xlabel = "",
        ylabel = "",
        title = "Zoom"
    )
    if template 
        savefig(folder_path*"/"*query*"_AF_CLUSTER_template/rmsd_scatter_"*query*"_U10.png")  # Sauvegarde du plot
    else 
        savefig(folder_path*"/"*query*"_AF_CLUSTER/rmsd_scatter_"*query*"_U10.png")  # Sauvegarde du plot
    end
        model_files = String[]
    # Trouver tous les dossiers qui commencent par "cluster"
    cluster_dirs = filter(isdir, glob("EX_U100*", folder_model))
    for cluster in cluster_dirs
        if template 
            models_path = cluster * "/af/predictions/" * split(basename(cluster), ".")[1] * "/models"  # 📂 Chemin des modèles
        else 
            models_path = joinpath(cluster, "models")
        end
        println(models_path)
        if isdir(models_path)   # Vérifier si le chemin existe
            matching_pdb = filter(isfile, glob("*.pdb", models_path))
            append!(model_files, matching_pdb)
        end
    end

    println(length(model_files))
    println(typeof(model_files))
    results=compare_model(model_files,holo_pdb,apo_pdb,folder_path,folder_model)
    if template 
        CSV.write(folder_path*"/"*query*"_AF_CLUSTER_template/rmsd_results_"*query*"_U100.csv", results)
    else 
        CSV.write(folder_path*"/"*query*"_AF_CLUSTER/rmsd_results_"*query*"_U100.csv", results)
    end
    #Visualisation of the result 
    
   scatter(
        results.RMSD_Apo,
        results.RMSD_Holo,
        xlabel = "Apo",
        ylabel = "Holo",
        title = "$query with AF-CLUSTER U100",
        marker = :circle,
        legend = false,
        xlims = (0, 30),
        ylims = (0, 30)
    )

    # Création de l'inset (zoom sur 0–7)
    scatter!(
        results.RMSD_Apo,
        results.RMSD_Holo,
        xlims = (0, 7),
        ylims = (0, 7),
        inset = (1, bbox(0.05, 0.05, 0.5, 0.5, :bottom, :right)),
        subplot = 2,
        marker = :circle,
        markersize = 3,  
        legend = false,
        xlabel = "",
        ylabel = "",
        title = "Zoom"
    )
    if template 
        savefig(folder_path*"/"*query*"_AF_CLUSTER_template/rmsd_scatter_"*query*"_U100.png")  # Sauvegarde du plot
    else 
        savefig(folder_path*"/"*query*"_AF_CLUSTER/rmsd_scatter_"*query*"_U100.png")  # Sauvegarde du plot
    end
    model_files = String[]
    cluster_labels = String[]  # Stocke le cluster d'origine pour chaque fichier
    cluster_dirs = filter(f -> occursin(r"^EX_(?!U)", basename(f)), glob("EX_*", folder_model))

    for cluster in cluster_dirs
        if template 
            models_path = cluster * "/af/predictions/" * split(basename(cluster), ".")[1] * "/models"  # 📂 Chemin des modèles
        else 
            models_path = joinpath(cluster, "models")
        end
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
    if template 
        CSV.write(folder_path * "/"*query*"_AF_CLUSTER_template/rmsd_results_"*query*".csv", results)
    else 
        CSV.write(folder_path * "/"*query*"_AF_CLUSTER/rmsd_results_"*query*".csv", results)
    end
    # Attribution d'une couleur unique par cluster
    unique_clusters = unique(cluster_labels)
    color_map = Dict(cluster => getindex(distinguishable_colors(length(unique_clusters)), i) for (i, cluster) in enumerate(unique_clusters))
    colors = [color_map[cluster] for cluster in cluster_labels]  # Associe la couleur à chaque point

    # Visualisation du résultat avec couleurs spécifiques
    
   scatter(
        results.RMSD_Apo,
        results.RMSD_Holo,
        xlabel = "Apo",
        ylabel = "Holo",
        title = "$query with AF-CLUSTER",
        marker = :circle,
        legend = false,
        xlims = (0, 30),
        ylims = (0, 30)
    )

    # Création de l'inset (zoom sur 0–7)
    scatter!(
        results.RMSD_Apo,
        results.RMSD_Holo,
        xlims = (0, 7),
        ylims = (0, 7),
        inset = (1, bbox(0.05, 0.05, 0.5, 0.5, :bottom, :right)),
        subplot = 2,
        marker = :circle,
        markersize = 3,  
        legend = false,
        xlabel = "",
        ylabel = "",
        title = "Zoom"
    )
    if template 
        savefig(folder_path * "/"*query*"_AF_CLUSTER_template/rmsd_scatter_"*query*".png")  # Sauvegarde du plot
    else 
        savefig(folder_path * "/"*query*"_AF_CLUSTER/rmsd_scatter_"*query*".png")  # Sauvegarde du plot
    end
end
################### MAIN #######################

folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
df_info = CSV.read(folder_path*"/info_dev_set.csv", DataFrame, delim=',')
println(df_info)
create_file=false
template=true

for row in eachrow(df_info)
    continuer =true
    query=row.PDB_apo
    println(query)
    ## AF-Cluster
    if template && create_file
        Alpha_conformer_result_path=joinpath(folder_path,query)
        cluster_dirs_AC = filter(isdir, glob("cluster*", Alpha_conformer_result_path))
        println(length(cluster_dirs_AC))
        if length(cluster_dirs_AC)!=0
            AF_cluster_template(query,create_file)
        else 
            println("We don't have template for $query")
        end
    elseif !template && create_file
        AF_cluster(query,create_file)
    else
        if template 
            folder_model=folder_path*query*"_AF_CLUSTER_template"
            if isdir(folder_model)
                folder_model=folder_path*query*"_AF_CLUSTER_template/output/"   
            else 
                println("We don't have the output of AF2 for $query")
                continuer=false
            end
        else 
            folder_model=folder_path*query*"_AF_CLUSTER/output/predictions/"   
        end
        if continuer
            filtered_rows = filter(row -> occursin(query, row.PDB_apo), df_info)
            println(filtered_rows)
            row = first(filtered_rows, 1)  # Prend la première ligne

            apo_pdb = string(row.PDB_apo[1], "_", row.CHAIN_apo[1], "_", row.INDEX_apo[1], ".pdb.gz")
            holo_pdb = string(row.PDB_holo[1], "_", row.CHAIN_holo[1], "_", row.INDEX_holo[1], ".pdb.gz")
        
            visualisation(folder_model,holo_pdb,apo_pdb,folder_path,query,template)
        end

    end
end
