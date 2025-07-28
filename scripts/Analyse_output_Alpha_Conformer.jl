#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analyse_output_Alpha_Conformer.jl.o%j.out
=#
import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")


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
function analyse_output(query::String,folder_path::String,folder_model::String)
    #Get the dev set file
    df_info = CSV.read(folder_path*"/info_test_set.csv", DataFrame, delim=',')
    filtered_rows = filter(row -> occursin(query, row.PDB_apo), df_info)
    row = first(filtered_rows, 1)  # Take the first row

    #Extract the apo and holo shape
    apo_pdb = string(row.PDB_apo[1], "_", row.CHAIN_apo[1], ".pdb")
    holo_pdb = string(row.PDB_holo[1], "_", row.CHAIN_holo[1], ".pdb")

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
    xlims=(0,25), ylims=(0, 25))
    #Save the plot
    savefig(folder_model*"/rmsd_scatter_"*query*".png")  # Sauvegarde du plot
end

"""
Analyse the output of AlphaConformer
Register the number of cluster, the number of file use , the lenght of the MSA...
"""
#Analyse the output directory of AlphaConformer 
function analyze_directory(base_dir::String,div::Bool)
    #Get all the cluster
    cluster_dirs = filter(isdir, sort(glob("cluster_*", base_dir)))
    #initialize the variabble
    n_clusters = length(cluster_dirs)
    template_counts = Int[]
    a3m_line_counts = Int[]
    template_files_set = Set{String}()
    all_seq_names = Set{String}() 

    #For each cluster
    for cluster_dir in cluster_dirs
        #Get the template
        template_dir = joinpath(cluster_dir, "templates")
        #Get the MSA
        a3m_files = glob("*.a3m", cluster_dir)
        
        # Count the number of file in template
        if isdir(template_dir)
            templates = readdir(template_dir)
            push!(template_counts, length(templates))
            foreach(f -> push!(template_files_set, f), templates)
        else
            push!(template_counts, 0)
        end

        # Count the number of row in MSA
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

    #Save the result
    return (
        n_clusters = n_clusters,
        avg_templates_per_cluster = avg_templates,
        avg_lines_per_a3m = avg_a3m_lines,
        template_files = template_files_set,
        unique_sequence_names = collect(all_seq_names)
    )
end

"""
Associates each template with the cluster number to which it belongs.
It returns a dictionary where each key is the name of a template file and the value is the cluster index.
"""
#Associate template to each cluster
function get_cluster_assignments(base_dir::String)
    cluster_dirs = filter(isdir, sort(glob("cluster_*", base_dir)))
    assignments = Dict{String, Int}()
    for (i, cluster_dir) in enumerate(cluster_dirs)
        template_dir = joinpath(cluster_dir, "templates")
        if isdir(template_dir)
            for f in readdir(template_dir)
                assignments[f] = i  # ou parse(Int, basename(cluster_dir)[end]) si cluster_X
            end
        end
    end
    return assignments
end

"""
Compares two template clusterings from two different folders.
It uses get_cluster_assignments to obtain the template assignments to each cluster in each folder, 
Then constructs two lists of cluster labels aligned with the set of templates in either folder.
It then calculates the Adjusted Rand Index (ARI)
Measure of similarity between two clusterings, and displays the result.
"""
#compares two template clusterings from two different folders.
function compare_clusterings(dir1::String, dir2::String)
    assign1 = get_cluster_assignments(dir1)
    assign2 = get_cluster_assignments(dir2)
    all_templates = union(keys(assign1), keys(assign2))
    labels1 = Int[]
    labels2 = Int[]
    for t in all_templates
        push!(labels1, get(assign1, t, -1))
        push!(labels2, get(assign2, t, -1))
    end
    ari = randindex(labels1, labels2)
    println("Adjusted Rand Index (ARI) between clusterings: $ari")
    return ari
end

"""
Take in input two directory of the output of AlphaConformers
Compare the result 
Use analyze_directory and compare_clusterings
Use it to compare first result and today result 
"""
#Compare output from two result of AlphaConformers
function compare_directories(dir1::String, dir2::String)
    stats1 = analyze_directory(dir1,true)
    stats2 = analyze_directory(dir2,false)

    println("=== Comparaison de $dir1 et $dir2 ===")
    println("Clusters: $(stats1.n_clusters) vs $(stats2.n_clusters)")
    println("Fichiers templates moyens: $(stats1.avg_templates_per_cluster) vs $(stats2.avg_templates_per_cluster)")
    println("Lignes moyennes .a3m: $(stats1.avg_lines_per_a3m) vs $(stats2.avg_lines_per_a3m)")

    files1_lower = Set(unique(filter(f -> occursin(r"\.", f), lowercase.(collect(stats1.template_files)))))
    files2_lower = Set(unique(filter(f -> occursin(r"\.", f), lowercase.(collect(stats2.template_files)))))

    only_in_dir1 = setdiff(files1_lower, files2_lower)
    only_in_dir2 = setdiff(files2_lower, files1_lower)

    println("\nFichiers uniquement dans $dir1/templates:")
    println(length(only_in_dir1))
    println(length(files1_lower))
    println("\nFichiers uniquement dans $dir2/templates:")
    println(length(only_in_dir2))
    println(length(files2_lower))

    # Overlap between sequence names (template names)
    seqs1 = Set(lowercase.(stats1.unique_sequence_names))
    seqs2 = Set(lowercase.(stats2.unique_sequence_names))
    intersection = intersect(seqs1, seqs2)
    union_set = union(seqs1, seqs2)
    iou = isempty(union_set) ? 0.0 : length(intersection) / length(union_set)

    println("\n=== Overlap des noms de séquence (IoU) ===")
    println("Nombre de séquences uniques dans $dir1 : $(length(seqs1))")
    println("Nombre de séquences uniques dans $dir2 : $(length(seqs2))")
    println("Intersection : $(length(intersection))")
    println("Union : $(length(union_set))")
    println("Intersection over Union (IoU) : $iou")

    compare_clusterings(dir1, dir2)
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
folder_path ="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"

#dev set that we want to look at
file_path="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/info_dev_set.csv"
info_pdb=DataFrames.DataFrame(CSV.File(file_path,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(info_pdb,20))
println(size(info_pdb))

#Files names that we want to look at 
file=["6akm","6r17","6uui","7lp1","9g21"] 
#OR file name 
query="6PWK"
################################################################################

# For multiple run 
for query in eachrow(info_pdb).PDB_apo
    query = string(query)  # Convert to string if necessary
    println("Analyse for query: ", query)
    #Folder output path from ALphaConformers
    folder_model = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/Devset_test/"*query*"_AlphaConformer_Hobohm"
    #Get the comparaison between every model and apo and hlo shape 
    analyse_output(query,folder_path,folder_model)
end
#=
#For one run 
#Get the comparaison between every model and apo and hlo shape 
analyse_output(query,folder_path,folder_model)

#Compare two AlphaConformer result 
########################## Information to fill #################################
#Path for the two folder
chemin_dossier_1 = joinpath(folder_path,"1AKZ_A/clusters")
chemin_dossier_2 = joinpath(folder_path,"1AKZ_Update_Hobohm_PDB")
###############################################################################
resultats = compare_directories(chemin_dossier_1, chemin_dossier_2)
@show "End"
########################################### END ##############################################################

=#