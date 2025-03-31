#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#PBS -l host=node48
#PBS -l walltime=900:00:00
#PBS -l mem=100gb
#PBS -l ncpus=40
#PBS -j oe
=#

#1AKZ
function AF_cluster()
    ### Create the folder to save the result
    folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
    cluster_dir=folder_path*"1AKZ_AF_CLUSTER_template"
    if !isdir(cluster_dir)
        mkdir(cluster_dir)
    end

    msa=joinpath(cluster_dir, "msa.a3m")
    if !isfile(msa)
        msa_origin=folder_path*"/1AKZ/afdb_up_results/msa.a3m"
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

    output_dir=cluster_dir*"/output"
    template=cluster_dir*"/template"
    COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
    af_command = `$COLABFOLD_PATH $msas $output_dir --use-templates 1 --msa-input --custom-template-path $template --overwrite-existing-results`
    @info "Running AlphaFold command: $af_command"
    run(af_command)
end

AF_cluster()
"""

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


## 1AKZ ###

#1AKZ,A,1,1SSP,E,1
folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
apo_pdb="1AKZ_A_1.pdb.gz"
holo_pdb="1SSP_E_1.pdb.gz"
folder_model=folder_path*"1AKZ_AF_CLUSTER/output/predictions/"

model_files = String[]
# Trouver tous les dossiers qui commencent par "cluster"
cluster_dirs = filter(isdir, glob("EX_U10*", folder_model))
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
CSV.write(folder_path*"/1AKZ_AF_CLUSTER/rmsd_results_1AKZ_U10.csv", results)

#Visualisation of the result 
scatter(results.RMSD_Holo, results.RMSD_Apo,
    xlabel="Holo", ylabel="Apo",
    title="Compare model with Apo and Holo form for 1AKZ with AF-CLUSTER U10",
    marker=:circle, legend=false,
    xlims=(0, 7), ylims=(0, 7))

savefig(folder_path*"/1AKZ_AF_CLUSTER/rmsd_scatter_1AKZ_U10.png")  # Sauvegarde du plot

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
CSV.write(folder_path*"/1AKZ_AF_CLUSTER/rmsd_results_1AKZ_U100.csv", results)

#Visualisation of the result 
scatter(results.RMSD_Holo, results.RMSD_Apo,
    xlabel="Holo", ylabel="Apo",
    title="Compare model with Apo and Holo form for 1AKZ with AF-CLUSTER U100",
    marker=:circle, legend=false,
    xlims=(0, 7), ylims=(0, 7))

savefig(folder_path*"/1AKZ_AF_CLUSTER/rmsd_scatter_1AKZ_U100.png")  # Sauvegarde du plot

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
CSV.write(folder_path * "/1AKZ_AF_CLUSTER/rmsd_results_1AKZ.csv", results)

# Attribution d'une couleur unique par cluster
unique_clusters = unique(cluster_labels)
color_map = Dict(cluster => getindex(distinguishable_colors(length(unique_clusters)), i) for (i, cluster) in enumerate(unique_clusters))
colors = [color_map[cluster] for cluster in cluster_labels]  # Associe la couleur à chaque point

# Visualisation du résultat avec couleurs spécifiques
scatter(results.RMSD_Holo, results.RMSD_Apo,
    xlabel="Holo", ylabel="Apo",
    title="Compare model with Apo and Holo form for 1AKZ with AF-CLUSTER",
    marker=:circle, legend=true,
    xlims=(0, 7), ylims=(0, 7),
    color=colors, label=cluster_labels)

savefig(folder_path * "/1AKZ_AF_CLUSTER/rmsd_scatter_1AKZ.png")  # Sauvegarde du plot
"""