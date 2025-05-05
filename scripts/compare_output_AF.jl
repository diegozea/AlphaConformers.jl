using DataFrames, CSV
using Statistics

#Our folder
path_main_dir="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"

#Get the pdb information
df_info_path=joinpath(path_main_dir,"Supplementary_Table_1_91_apo_holo_pairs.csv")
df_info=DataFrames.DataFrame(CSV.File(df_info_path,
        comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
#Create the DF
df_comparaison = DataFrame(Folder=String[],File = String[], rmsd_apo_holo = Float64[], min_apo=Float64[], mean_apo =Float64[], var_apo=Float64[],min_holo=Float64[], mean_holo = Float64[], var_holo =Float64[])

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
            row_match = filter(row -> occursin(parts[1], row.apo_id), df_info)
            rmsd_apo_holo=row_match.rmsd_apo_holo[1]

            #Get the file 
            df_file=DataFrames.DataFrame(CSV.File(csv_file,
            comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
            push!(df_comparaison,(folder_name,file_name,rmsd_apo_holo,minimum(df_file.RMSD_Apo),mean(df_file.RMSD_Apo),var(df_file.RMSD_Apo),minimum(df_file.RMSD_Holo),mean(df_file.RMSD_Holo),var(df_file.RMSD_Holo)))
        catch e
            println("didn't work for $folder_name",e)
            continue
        end
        end
end
CSV.write("Compare_output_AF.csv", df_comparaison)
println(first(df_comparaison,20))