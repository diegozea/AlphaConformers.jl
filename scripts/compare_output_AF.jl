# Load necessary packages
using DataFrames, CSV
using Statistics
using Glob
using MIToS
using AlphaConformers

############################################################ Functions #############################################################"
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

##################################################### MAIN #############################################################
"""
Code to save the performance of multiple test 

Can be use on AlphaConformer or other like AF-cluster and BioEmu 
Need juste the csv file that register the rmsd between apo and holo 
--> directly created with BIOEMU_development_test.jl and AF_cluster_development_test.jl
--> run Analyse_output_Alpha_Conformer.jl for AlphaConformer 

Input :
- path_main_dir: path to the directory that have all the result 
--> will go in every folder search for the csv file 
Output : 
CSV file "Compare_output_AF_test_set.csv" thta save the best performance for each test 

Create the plot with compare_output_AF_visualization.jl
"""
########################## Information to fill #################################
#Path to the folder with all the result 
path_main_dir="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/test_set/"

#Get the pdb information
df_info_path=joinpath(path_main_dir,"info_test_set.csv")

################################################################################

df_info=DataFrames.DataFrame(CSV.File(df_info_path,
        comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
#Create the DF
df_comparaison = DataFrame(Folder=String[],File = String[], rmsd_apo_holo = Float64[], Nombre_cluster = Float64[], Mean_MSA= Float64[], Sequence_in_msa = Int128[], min_apo=Float64[], mean_apo =Float64[], var_apo=Float64[],min_holo=Float64[], mean_holo = Float64[], var_holo =Float64[])
 
#get all the folders 
folders = filter(isdir, readdir(path_main_dir, join=true))

#For each folder
for folder in folders
    folder_name=basename(folder)
    csv_files = filter(f -> isfile(f) && endswith(f, ".csv"), readdir(folder, join=true))
    for csv_file in csv_files
        try 
            println(basename(csv_file))
            #Take the file name 
            filename_no_ext = split(basename(csv_file), ".")[1]
            println(filename_no_ext)
            file_name = join(split(filename_no_ext, "_")[3:end], "_")
            println(file_name)

            #Get the rmsd 
            parts = split(file_name, "_")
            row_match = filter(row -> occursin(parts[1], row.PDB_apo), df_info)
            rmsd_apo_holo=row_match.RMSD[1]

            #Get the file 
            result_csv_AF=DataFrames.DataFrame(CSV.File(csv_file,
            comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
            result_directory=analyze_directory(folder,false)
            result_directory_seqs = Set(lowercase.(result_directory.unique_sequence_names))
            push!(df_comparaison,(folder_name,file_name,rmsd_apo_holo,result_directory.n_clusters,result_directory.avg_lines_per_a3m,length(result_directory_seqs),minimum(result_csv_AF.RMSD_Apo),mean(result_csv_AF.RMSD_Apo),std(result_csv_AF.RMSD_Apo),minimum(result_csv_AF.RMSD_Holo),mean(result_csv_AF.RMSD_Holo),std(result_csv_AF.RMSD_Holo)))
            
        catch e
            println("didn't work for $folder_name",e)
            continue
        end
    end
end

CSV.write("Compare_output_AF_test_set.csv", df_comparaison)
@show first(df_comparaison,20)

@show "END"

####################################################################### END #########################################################