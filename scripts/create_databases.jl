#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=20
#SBATCH --output=create_databases-%j.out
=#

import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")
using Distributed
addprocs(19)  # Ajoute 4 processus parallèles

# Load necessary packages on all workers
@everywhere using DataFrames, CSV
@everywhere using Dates
@everywhere using MIToS.PDB
@everywhere import MIToS
@everywhere using AlphaConformers
@everywhere using BioStructures
@everywhere using DataStructures
@everywhere using Statistics
@everywhere using Glob
@everywhere using Foldseek_jll


############################################### Fonctions ######################################################

# Check that apo and holo are in different cluster
@everywhere function filter_pdb(df_uniprot,pdb_folder::String)
    # Get all the apo and holo conformation
    with_ligands = filter(row -> !ismissing(row.LIGANDS), df_uniprot)
    without_ligands = filter(row -> ismissing(row.LIGANDS), df_uniprot)
    result = 0

    #Check every cluster
    for cutoff in [0.5, 0.75, 1.0, 1.25, 1.5,2]
        cluster_col = Symbol("Cluster_$(cutoff)")

        for file_with_ligands in eachrow(with_ligands) # For the holo conformation
            for file_without_ligands in eachrow(without_ligands) # For the apo conformation
                cluster_with = file_with_ligands[cluster_col]
                cluster_without = file_without_ligands[cluster_col]
                if !ismissing(cluster_with) && !ismissing(cluster_without) && cluster_with != cluster_without
                    # If the files are in different clusters, take the smaller value
                    if  cutoff>result # We want the maximum cutoff
                        result = cutoff
                    end 
                end
            end
        end
    end

    # If we have a result
    if result != 0
        #check that they have more than 1A of differences 
        if result >=1.0

            cutoff =1.0
            cluster_col = Symbol("Cluster_$(cutoff)")
           # Extract the apo and holo conformations for cutoff 1.0
            clusters_with = Set(row[cluster_col] for row in eachrow(with_ligands) if !ismissing(row[cluster_col]))
            clusters_without = Set(row[cluster_col] for row in eachrow(without_ligands) if !ismissing(row[cluster_col]))

            # Keep only the clusters values that are not in both apo and holo conformations
            if isempty(intersect(clusters_with, clusters_without))
                df_final = vcat(with_ligands, without_ligands)
            else
                df_final = DataFrame()  # empty DataFrame if all clusters match
            end
            # If we have a final DataFrame, check that we still have apo
            if !isempty(df_final)
                #check that we still have apo 
                is_apo = any(ismissing, df_final.LIGANDS)
                if is_apo
                    return df_final
                else 
                    return nothing
                end
            else 
                return nothing
            end
        else
            return nothing
        end
    else 
        return nothing
    end
    
end


"""
This function takes in input a DataFrame with the columns UNIPROT, PDB and CHAIN
It returns a DataFrame with the columns UNIPROT, PDB, CHAIN
It uses the function filter_pdb to filter the PDB files and keep only those that have apo and holo conformations in different clusters
It uses pmap to parallelize the processing of the DataFrame
It processes the DataFrame in chunks of 10 groups at a time
"""
#Create a database with the PDB files that have apo and holo conformations in different clusters
@everywhere function create_database(df_final::DataFrame,pdb_folder:: String)
    grouped_uniprots = groupby(df_final, :UNIPROT)
    chunks = Iterators.partition(grouped_uniprots, 10)  # 10 groups at a time

    # process each group
    results = pmap(chunk -> begin
        tmp = Vector{DataFrame}()
        for group in chunk
            res = filter_pdb(group,pdb_folder) #check if the apo and holo conformations are in different clusters
            if res !== nothing
                push!(tmp, res) # Add the result to the temporary vector
            end
        end
        tmp
    end, collect(chunks))

    return vcat(Iterators.flatten(results)...) # Combine all results into a single DataFrame
end


"""
Use the function getpdbdescription to have the details about the pdb file 

Take in input a PDB code in lowercase and output the information in a vecteur 
We will have the code pdb, the release date, the method of extraction and the resolution if we have the information 
(can be missing if the method is NMR)
"""
# Function to get the header information of a PDB file
@everywhere function get_header_pdb(code_pdb::String)
    header_info=getpdbdescription(code_pdb) # Returns a DICT

    #Extract the information
    ## Retrieve the publication date
    date_info_p = get(header_info, "rcsb_accession_info", Dict())
    date_info_d = get(date_info_p, "initial_release_date", nothing)
    if date_info_d !== nothing
        date_info_t = split(date_info_d, "T")
        date_info = date_info_t[1]
        return date_info
    else
        return nothing
    end   
end

"""
It takes in input a DataFrame with the columns UNIPROT, PDB and CHAIN
It returns a DataFrame with the columns UNIPROT, PDB, CHAIN and NB_SIMILAR_PROT     
It also returns a DataFrame with the columns UNIPROT, PDB and CHAIN for the proteins that have no similar structure
It uses Foldseek to find similar structures
"""
# Function to find similar PDB files using Foldseek
@everywhere function foldseek_similar_pdb(df_final::DataFrame, FOLDSEEK_DB::String, pdb_folder::String, details_df::DataFrame)
    ##################
    # Function that processes a single line
    @everywhere function process_row(row,FOLDSEEK_DB,pdb_folder)
        # Extract the UNIPROT, PDB and CHAIN from the row
        uniprot = row.UNIPROT
        pdb_file = row.PDB
        chain = row.CHAIN
        database_row = Vector{Any}() 
        @show "Processing $uniprot - $pdb_file - $chain"
        # Check if PDB file exists
        pdbfilepath = joinpath(pdb_folder, uppercase(pdb_file) * ".pdb")
        @show pdbfilepath
        # If the PDB file exists
        if isfile(pdbfilepath)
            try
                #launch Foldseek to find similar structures
                @info "Run Foldseek for $pdb_file"
                foldseek_df = mktempdir() do tmp_folder
                    result_file = joinpath(tmp_folder, "foldseek_result.m8")
                    run(`$(Foldseek_jll.foldseek()) easy-search $pdbfilepath $FOLDSEEK_DB $result_file $tmp_folder --format-output query,target,evalue,alntmscore`)
                    if isfile(result_file)
                        foldseek_df = CSV.read(result_file, DataFrame, delim='\t', header=[
                            "query", "target", "evalue","tmscore"
                        ])
                        return foldseek_df
                    else
                        return DataFrame()
                    end
                end
                if isempty(foldseek_df)
                    database_row= (UNIPROT=uniprot, PDB=pdb_file, CHAIN=chain, NB_SIMILAR_PROT=0)
                    return database_row # Return the row with 0 similar structures
                end
                @show first(foldseek_df)
                # Filter results to keep those with evalue < 1e-5 and tmscore > 0.5
                df_result_evalue = filter(r -> r.evalue < 1e-5 && r.tmscore>0.5, foldseek_df)

                count = 0
                rows_to_keep = Vector{Int}()
                # Check results output and retrieve publication date
                for (i, r) in enumerate(eachrow(df_result_evalue))
                    # Extract target and chain from the result
                    target = String(split(r.target, ".")[1])
                    parts = split(r.target, "_")
                    release_date = missing
                    #Get the publication date
                    if length(parts) >= 2
                        chain2 = parts[2]
                        # Get the rows with the same PDB and chain in details_df
                        # We use lowercase to match the PDB code
                        # and chain2 to match the chain
                        matching_rows = filter(rr -> rr.PDB == lowercase(target) && rr.CHAIN == chain2, eachrow(details_df)) 
                        if !isempty(matching_rows) # If we have matching rows
                            release_date = matching_rows[1][:DATE_RELEASE] # Extract the release date
                        else
                            continue # If no matching rows, skip this iteration
                            # count as no similar structure found because release after the creation of the database (2025-03-01)
                        end
                    else
                        continue    
                    end

                    # Check if the release date is before 2018-04-30
                    if !ismissing(release_date) && Date(release_date) < Date("2018-04-30")
                        push!(rows_to_keep, i)
                        count += 1 #count as a similar structure found
                    end
                end
                #Add information to the database row
                database_row= (UNIPROT=uniprot, PDB=pdb_file, CHAIN=chain, NB_SIMILAR_PROT=count) 
            catch e
                @warn "Erreur sur $pdb_file : $e"
            end
        else
            @warn "The PDB file of $pdb_file doesn't exist"
        end
        return database_row # return the database row
    end
    ###########################
    
    # Initialize the final DataFrame to store results
    database_final = DataFrame(UNIPROT=String[], PDB=String[], CHAIN=String[], NB_SIMILAR_PROT=Int64[])
    # Lire le fichier intermédiaire si il existe
    intermediate_path = "save_intermediate_result.csv"
    database_final = isfile(intermediate_path) ? DataFrame(CSV.File(intermediate_path)) : DataFrame(UNIPROT=String[], PDB=String[], CHAIN=String[], NB_SIMILAR_PROT=Int[])
    @show first(database_final,10)
    # Filtrer df_final pour ne garder que les PDB/CHAIN non déjà traités
    df_to_process = filter(row -> !(row.PDB in database_final.PDB && row.CHAIN in database_final.CHAIN), df_final)
    @show first(df_to_process, 20)
    @show size(df_to_process)
    # Créer les chunks à partir de df_to_process
    chunks = Iterators.partition(eachrow(df_to_process), 10)

    #chunks = Iterators.partition(eachrow(df_final), 10)  # Process in chunks of 10 rows 
    @info "Processing rows in chunks"
    # Process each chunk in parallel
    for chunk in chunks
        @info "Processing chunk"
        # Process each row in the chunk in parallel
        database_row = pmap(row -> process_row(row,FOLDSEEK_DB,pdb_folder), collect(chunk))
        # Combine the results into a single DataFrame
        @info "Combining results"
        df_completed = vcat(skipmissing(database_row)...)
        append!(database_final, df_completed)
        CSV.write("save_intermediate_result.csv", database_final, append=false) # Save results to CSV file
        GC.gc()  # Run garbage collection to free memory
    end
    
        
    # Process the results to get the final DataFrame
    @info "Results processed"
    df_merged = innerjoin(df_final, database_final, on=[:UNIPROT, :PDB, :CHAIN])
    return df_merged # return the final DataFrame with NB_SIMILAR_PROT
end


########################################## MAIN #######################################################"
"""
Main function to create the database
It reads the PDB files, filters them, and finds similar structures using Foldseek
It saves the final DataFrame to a CSV file

Input:
- file_path_df_final: path to the CSV file with PDB information output by mapping_Uniprot_PDB.jl
- pdb_folder: String with the path to the database of PDB 
- FOLDSEEK_DB: String with the path to the Foldseek database
- filtering: Boolean to filter the DataFrame or not to keep only the PDB files with apo and holo conformations in different clusters
- save: Boolean to save all the DataFrame to a CSV file
- details_path: String with the path to the CSV file with the all details of all the PDB files
- filtering2: Boolean to filter the DataFrame or not to keep only the PDB files with no similar structures and with apo and holo conformations in different clusters
Output:
- df_merged: DataFrame with PDB information and no similar structures found before 2018-04-30 

Take aroound 5 hours for each step
Advise to put save to true to save the DataFrame at each step
Be careful of the memory usage, it can take a lot of memory
"""
########################## Information to fill #################################

# Define the path to the Foldseek database
const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

# Define the path to the PDB folder
const pdb_folder= "/alpha/database/pdb/pdb_files"

# Get the file path with the result from mapping_Uniprot_PDB.jl
file_path_df_final="pdb_information_details_final_mutation_cluster_reformatted.csv" 

# Filter the DataFrame to keep only the PDB files with apo and holo conformations in different clusters
filtering = false

# Save each file for each step
save = true

# Get the details of all the PDB files
details_path="pdb_information_details.csv" 

# Filter the DataFrame to keep only the PDB files with no similar structures and with apo and holo conformations in different clusters
filtering2 = true
##################################################################################

@info "Start"
@show "Date de debut ", Dates.format(now(), "HH:MM:SS")

# Read the DataFrame from the CSV file
df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(df_final,20))
println(size(df_final))

if filtering 
    # Filter the DataFrame to keep only the PDB files with apo and holo conformations in different clusters
    df_final=create_database(df_final,pdb_folder)
    #if save the DataFrame is save in a CSV file
    if save 
        CSV.write("pdb_information_details_final_mutation_cluster_reformatted_filter.csv", filter(!isnothing, df_final)) 
    end
end
println(first(df_final,20))
println(size(df_final))


if save 
    # If we save the DataFrame, we read it again to have the DataFrame
    file_path_df_final="pdb_information_details_final_mutation_cluster_reformatted_filter.csv" # File path
    df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
    comment="#", missingstring=["", "None"]))
    @show size(df_final)
end

# Read the DataFrame from the CSV file with the details of the PDB files
details_df = DataFrame(CSV.File(details_path, comment="#", missingstring=["", "None"]))
    
if df_final !== nothing 
    #Get the number of similar structures using Foldseek
    pdbfilepath = joinpath(pdb_folder, "1AKZ.pdb")
    @show isfile(pdbfilepath)
    df_merged =foldseek_similar_pdb(df_final,FOLDSEEK_DB,pdb_folder,details_df)
    println(first(df_merged,20))
    println(size(df_merged))
    if save 
        CSV.write("pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek.csv",df_merged)
    end
end

if save 
    # If we save the DataFrame, we read it again to have the DataFrame
    file_path_df_merged="pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek.csv" # File path
    df_merged=DataFrames.DataFrame(CSV.File(file_path_df_merged,
    comment="#", missingstring=["", "None"]))
    @show size(df_merged)
end

if filtering2 
    # Filter the DataFrame to keep only the PDB files with no similar structures
    df_merged = filter(row -> row.NB_SIMILAR_PROT == 0, df_merged)
    # If we have no similar structures, we keep the PDB files with apo and holo conformations in different clusters
    df_final=create_database(df_merged,pdb_folder)

    println(first(df_final,20))
    println(size(df_final))

    #save the final DataFrame to a CSV file
    CSV.write("pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek_final.csv", filter(!isnothing, df_final)) 
end

@show "Date de fin ", Dates.format(now(), "HH:MM:SS")
@info "End"
############################################ END #######################################################
