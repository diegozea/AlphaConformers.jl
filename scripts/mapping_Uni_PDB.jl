#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l ncpus=13
#PBS -l mem=16g
#PBS -l host=node48
#PBS -l walltime=600:00:00
#PBS -j oe
=#

import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")
using Distributed
addprocs(12)  # Ajoute 12 processus parallèles

@everywhere using AlphaConformers
@everywhere using BioStructures
@everywhere using DataFrames
@everywhere using MIToS.PDB
@everywhere using MIToS.SIFTS
@everywhere using MIToS.MSA
@everywhere using JSON3
@everywhere using Statistics
@everywhere import CSV
@everywhere import MIToS
@everywhere using OrderedCollections
@everywhere using DataStructures
@everywhere using Plots
@everywhere import BioAlignments
@everywhere using Clustering 
@everywhere using Statistics
@everywhere using LinearAlgebra
@everywhere using PairwiseListMatrices
@everywhere using StatsBase
@everywhere using Dates
@everywhere using Profile


#Link the 3 csv together
"""
Create one file with all of the information form PFAM CATH and UNIPROT from SIFT databases 

Take in input the three file in a DataFrame format and output one DataFrame
"""
function join_information(sift_pfam_mapping::DataFrame,sift_cath_mapping::DataFrame,sift_uniprot_mapping::DataFrame)

    # 1 Extract RES_BEG and RES_END from sift_uniprot_mapping
    uniprot_res = select(sift_uniprot_mapping, [:PDB, :CHAIN, :RES_BEG, :RES_END, :SP_PRIMARY])
    # 2 Count the number of PFAM_IDs per (PDB, CHAIN, SP_PRIMARY)
    pfam_counts = combine(groupby(sift_pfam_mapping, [:PDB, :CHAIN, :SP_PRIMARY])) do subdf
        (PFAM_NB = size(subdf, 1),)
    end
    # 3 Count the number of CATH_IDs per (PDB, CHAIN, SP_PRIMARY)
    cath_counts = combine(groupby(sift_cath_mapping, [:PDB, :CHAIN, :SP_PRIMARY])) do subdf
        (CATH_NB = size(subdf, 1),)
    end
    # 4 Join the two tables to create the final DataFrame
    sift_join_file = outerjoin(uniprot_res, pfam_counts, on=[:PDB, :CHAIN, :SP_PRIMARY])
    sift_join_file = outerjoin(sift_join_file, cath_counts, on=[:PDB, :CHAIN, :SP_PRIMARY],)
    # 5 Replace `missing` values ​​with 0 (in case a PDB has no PFAM or CATH)
    sift_join_file.PFAM_NB .= coalesce.(sift_join_file.PFAM_NB, 0)
    sift_join_file.CATH_NB .= coalesce.(sift_join_file.CATH_NB, 0)
    # 6 Sort by SP_PRIMARY
    sift_join_file = sort(sift_join_file, :SP_PRIMARY)

    return sift_join_file 
end 

#Retrieve useful header information
"""
Use the function getpdbdescription to have the details about the pdb file 

Take in input a PDB code in lowercase and output the information in a vecteur 
We will have the code pdb, the release date, the method of extraction and the resolution if we have the information 
(can be missing if the method is NMR)
"""
function get_header_pdb(code_pdb::String)
    header_info=getpdbdescription(code_pdb) # Returns a DICT
    #Extract the information
    ## Retrieve the publication date
    date_info_p = get(header_info, "rcsb_accession_info", Dict())
    date_info_d = get(date_info_p, "initial_release_date", nothing)
    if date_info_d !== nothing
        date_info_t = split(date_info_d, "T")
        date_info = date_info_t[1]
    else
        return nothing
    end
    ## Retrieve method
    exp_info = get(header_info, "rcsb_entry_info", Dict())
    method_info = get(exp_info, "experimental_method", nothing)
    ## Retrieve the resolution
    res_info_p = get(exp_info, "resolution_combined", nothing)
    #Prepare the output
    if res_info_p !== nothing && !isempty(res_info_p)
        res_info = res_info_p[1]
        return (code_pdb, date_info, method_info, res_info)
    else 
        return (code_pdb, date_info, method_info,missing)
    end
end

#Save information in a DF
"""
From the Dataframe with all the information of Pfam CATH and Uniprot we had the information form the header 
We use the function get_header_pdb

Take in input the Dataframe from join_information and output a Dataframe with the header information added 
"""
function get_pdb_information(sift_join_file::DataFrame)
    #Create empty DF
    df_pdb_reso_resi = DataFrame(PDB = String[], CHAIN = String[], RES_BEG=Union{Int64,Missing}[], RES_END=Union{Int64,Missing}[],PFAM_NB=Int64[],CATH_NB=Int64[], RESOLUTION = Union{Float64, Missing}[], METHOD = String[],DATE_RELEASE=SubString{String}[],UNIPROT=String[]) # Take into account if missing resolution   
    #For each pdb file   
    for row in eachrow(sift_join_file)  
        # Retrieve associated uniprot code
        uni_acc=row.SP_PRIMARY
        #Retrieve header information 
        code_pdb, date_info, method_info, res_info = get_header_pdb(String(row.PDB))
        #Fill out the DF
        push!(df_pdb_reso_resi,(code_pdb,String(row.CHAIN),row.RES_BEG,row.RES_END,row.PFAM_NB,row.CATH_NB,res_info,method_info,date_info,String(uni_acc))) 
    end 
    return df_pdb_reso_res
end

#Recover the ligands
"""
Get the type of ligand that can bound with each pdb, if no ligand we put missing

Take in input the DataFrame from BioLip and the DataFrame created with get_pdb_information and add to it in output 
the ligand (We are looking the type not the quantities) 
"""
function get_ligand_information(df_biolip_5_first_columns::DataFrame,df_pdb_reso_resi::DataFrame)
    # Get a line for each pdb --> group the ligands
    df_ligands = combine(groupby(df_biolip_5_first_columns, [:PDB, :CHAIN]), 
                     :LIGAND => (x -> join(unique(x), ", ")) => :LIGANDS)
    # Merge with the main DF
    df_final = leftjoin(df_pdb_reso_resi, df_ligands, on=[:PDB, :CHAIN])
    # Replace missing values ​​with "missing"
    df_final.LIGANDS .= coalesce.(df_final.LIGANDS, missing)
    return df_final
end 

"""
Do the correspondance between Uniprot and PDB index. Use SIFTS mapping so we need to download the .xml file.
We handle the error if the file couldn't be install 

Take in input a row form the DataFrame get_ligand_information and the path to a folder to save temporary file
output a dictionnaire with the mapping Uniprot => PDB
"""
#Retrieve the link between Uniprot and PDB index
@everywhere function get_sift_mapping(row::DataFrameRow,sift_folder::String)
    pdb=row.PDB
    #Path where the xml file is saved 
    siftsfile=joinpath(sift_folder,pdb*".xml.gz") #Use the temporary folder
    if !isfile(siftsfile) #Check that the file is not already downloaded
        @warn "❌ Error Sift file of $pdb is not download "
        return nothing
    end
    #Do the mapping 
    siftsmap = siftsmapping(  # Returns a Dictionary with Uniprot coordinate => PDB coordinate
        siftsfile,
        dbUniProt,
        String(row.UNIPROT),
        dbPDB,
        String(pdb), # SIFTS stores PDB identifiers in lowercase
        chain = String(row.CHAIN),
        missings = false,
    ) # Residues without coordinates aren't used in the mapping
    return siftsmap #Output the mapping UNIPROT => PDB
end


"""
Make the link between the Uniprot index and the residues. Use the fonction get_sift_mapping and MIToS.PDB.downloadpdb

Take in input one row form the DataFrame get_ligand_information and the path to a folder to save temporary file
output a dictionnaire with the mapping Uniprot => Residues
"""
#Relating Uniprots indices to residuals
@everywhere function get_uniprot_mapping_residues(row::DataFrameRow,sift_folder::String,pdb_folder::String)
    # Retrieve the mapping between Uniprot => PDB
    mapping=get_sift_mapping(row,sift_folder) #OrderedDict Uniprot => PDB
    if mapping === nothing # If the xml file couldn't be download
        return nothing # don't do the mapping 
    end
    # Recover PDB residues
    pdbfile=joinpath(pdb_folder,uppercase(row.PDB)*".pdb" ) # Use the temporary file
    pdb=row.PDB
    if !isfile(pdbfile) #Check that the file is not already downloaded
        @warn "❌ Error the PDB file $pdb is not download"
        return nothing  
    end
    #Create the mapping
    residues_1ivo = read(pdbfile, PDBFile)
    chain_residues_pdb = MIToS.PDB.residuesdict(residues_1ivo; model="1", chain=String(row.CHAIN), group="ATOM") #Returns dictionary with position => residues details 
    #For all residues find the Uniprot index
    residues_uniprot = OrderedDict()
    for id_uni in keys(mapping)
        try 
            id_pdb=mapping[string(id_uni)] #Must be a string to compare them 
            residue=chain_residues_pdb[id_pdb]
            residues_uniprot[string(id_uni)]=residue # Output Ordred Dict Uniprot => Residue 
        catch e 
            pdb=row.PDB
            @warn "❌ Key $id_uni not found for $pdb" #If some residue aren't represent in the pdb we ignore them 
        end 
    end
    return residues_uniprot
end


"""
Compare the residues form PDB and Uniprot to identify the number of mutation and the missing residues

Take in input the DataFrame from get_ligand_information and the path to a folder to save temporary file
Output the Dataframe with two new colonnes MUTATION and MISSING_RESIDUES
"""
#Register the number of mutation and missing residues
@everywhere function  count_mutation_pdb_uni(df_uniprot::DataFrame,sift_folder::String)
    #Add column with mutation
    df_uniprot.MUTATION .= 0 
    #Add column with missing residue
    df_uniprot.MISSING_RESIDUES .= 0
    #For each pdb file
    for row in eachrow(df_uniprot)
        mut=0
        pdb=row.PDB
        @info "Uniprot " row.UNIPROT        
        #Recover the SIFT file
        siftsfile=joinpath(sift_folder,pdb*".xml.gz") #Use the temporary folder
        if !isfile(siftsfile)  # If the xml couldn't be download  
            @warn "❌ Error sift file $pdb is not download"
            return df_uniprot #Don't count the mutation and missing residues  
        end
        #Do the comparaison 
        residue_data = read(siftsfile, SIFTSXML)
        #For each residues
        for i in 1:length(residue_data)
            #get the information
            residues=read(siftsfile, SIFTSXML)[i]
            res_uni=get(residues, dbUniProt, :name, "missing")
            id_uni=get(residues, dbUniProt, :id, "missing")
            res_pdb=get(residues, dbPDB, :name, "")
            chain_pdb=get(residues, dbPDB, :chain, "") #Check that we are on the right channel and good Uniprot
            #If not in uniprot we ignore
            if row.CHAIN==chain_pdb && row.UNIPROT==id_uni
                res_pdb_n=three2residue(res_pdb) # Returns a MIToS.MSA.Residue
                if res_uni != string(res_pdb_n) # If not same residue then mutation
                    mut=mut+1
                end
            end  
        end
        #Add value in df
        row.MUTATION = mut
        #Count the number of missing residues
        sifts_1ivo = read(siftsfile, SIFTSXML, chain = String(row.CHAIN))
        mis=count([res.missing for res in sifts_1ivo]) # Returns a boolean vector so counts the number of 1s
        row.MISSING_RESIDUES=mis
    end
    return df_uniprot
end


"""
Compare for each uniprot accession number the pdb file associate with structure alignment 

Take in input the DataFrame from count_mutation_pdb_uni and the path to a folder to save temporary file
Output a Vector with the RMSD between each PDB --> its the half of the Correlation matrix
"""
#Compute the RMSD between each pdb 
@everywhere function calculate_RMSD_uniprot(df_uniprot::DataFrame,sift_folder::String,pdb_folder::String)
    files = df_uniprot.PDB .* "_" .* df_uniprot.CHAIN  # List of files
    n = length(files)  # Number of files
    rmsd_list = []
    coverage_list=[]
    for i in 1:n
        for j in (i+1):n  # Avoid duplicates
            #Get the mapping Uniprot => Residues for both file
            dict1=get_uniprot_mapping_residues(df_uniprot[i,:],sift_folder,pdb_folder) 
            dict2 = get_uniprot_mapping_residues(df_uniprot[j,:],sift_folder,pdb_folder)
            if dict1 === nothing  || dict2 === nothing # Check if the mapping was possible
                push!(coverage_list,0.0)
                rmsd = missing # we will not annalyse those pdb 
            else
                common_keys = Set(keys(dict1)) ∩ Set(keys(dict2)) # get the common key between the 2 pdb 
                #Calcul the pourcentage of coverage 
                pourcentage_coverage = length(common_keys) / max(length(dict1), length(dict2))
                push!(coverage_list,pourcentage_coverage)
                if pourcentage_coverage==0.0 # If no coverage 
                    rmsd=missing # we will not annalyse those pdb 
                else 

                    positions = SortedSet(parse.(Int,common_keys))  # Retrieves common positions (same key in all files)
                    # Extracting residues for each file as an ordered list
                    residues1 = [dict1[string(pos)] for pos in positions] 
                    residues2 = [dict2[string(pos)] for pos in positions]
                    try  
                        # Superposition and calculation of RMSD
                        _, _,rmsd = superimpose(residues1, residues2)
                    catch e 
                        @show common_keys
                        @warn "❌ Error while superimposing structures "
                        rmsd = missing # we will not annalyse those pdb 
                    end 
                end
            end
            push!(rmsd_list, rmsd)  # same the rmsd
        end
    end
    # Create the matrix of coverage
    plm=PairwiseListMatrix(coverage_list,false,1.0)
    row_means = [mean(row) for row in eachrow(plm)]
    keep_indices = findall(x -> x > 0.8, row_means)
    return rmsd_list,files,keep_indices # Output the file with all the rmsd and the file name
end


"""
Create the cluster for different CutOff : 1.5 2 and 4

Take in input the Vector with the RMSD from calculate_RMSD_uniprot
Output a Dataframe with the cluster number for each CutOff
"""
#Create the cluster for each uniprot accession number
@everywhere function cah(rmsd_list,files,keep_indices)
    #Create the symmetric matrix
    plm=PairwiseListMatrix(rmsd_list,false,0.0)
    #nplm=setlabels(plm,files)
    filtered_mat = plm[keep_indices, keep_indices] #To extract the file without the right pourcentage of coverage 
    # Delete rows and columns that contain `missing` --> where we couldn't analyse them
    rows_with_missing = findall(row -> any(ismissing, row), eachrow(filtered_mat))
    df_rmsd_2 = filtered_mat[setdiff(1:end, rows_with_missing), setdiff(1:end, rows_with_missing)]

    mat = Matrix{Float64}(df_rmsd_2) # Conversion to Matrix
    clustering_result = hclust(mat;linkage=:average)#Do the clustering

    #Cut the clusters according to different cutoffs
    cutoffs = [0.5, 0.75, 1, 1.25, 1.5,2] # List of different cutoffs to test

    #Creation of the output DF
    pdb_values = [split(file, "_")[1] for file in files]
    chain_values = [split(file, "_")[2] for file in files]
    final_df = DataFrame("PDB" => pdb_values,"CHAIN"=>chain_values)

    for cutoff in cutoffs  # Cut clusters at different cutoff levels
        cluster_assignments = cutree(clustering_result, h=cutoff)
        #Save the result
        df_clusters = DataFrame()
        file_i = files[setdiff(1:end, rows_with_missing)]
        df_clusters.PDB = [split(f, "_")[1] for f in file_i]  
        df_clusters.CHAIN = [split(f, "_")[2] for f in file_i] 

        cluster_col = Symbol("Cluster_$cutoff")
        if length(cluster_assignments) == length(file_i)
            df_clusters[!, cluster_col] = cluster_assignments
        else
            @info "⚠️ Cluster length mismatch, filling with `missing` for cutoff = $cutoff"
            @show length(cluster_assignments)
            @show cluster_assignments
            @show length(file_i)
            @show file_i
            @show rows_with_missing
            df_clusters[!, cluster_col] = Vector{Union{Missing, Int}}(missing, length(file_i))
        end
        #df_clusters[!, Symbol("Cluster_$cutoff")] = cluster_assignments  # Add the cluster to the DF

        # Do a `left join` to include PDBs that could not be parsed
        final_df = leftjoin(final_df, df_clusters, on=[:PDB,:CHAIN])
    end
    return final_df
end
#Filter the pdb to keep the relevant one 
@everywhere function filter_pdb(df_uniprot)
    #Check that we have apo and holo conformation 
    is_apo = is_apo = any(ismissing, df_uniprot.LIGANDS) && any(!ismissing, df_uniprot.LIGANDS)
    if is_apo 
        df_uniprot_res = filter(row -> ismissing(row.RESOLUTION) || row.RESOLUTION < 3, df_uniprot) #Pour recuperer fichier avec bonne résolution
        df_uniprot_res_date = filter(row -> row.DATE_RELEASE >= Date("2018-01-01", "yyyy-mm-dd"), df_uniprot_res) #get pdb release after the training of AF2
        println(df_uniprot_res_date)
        if !isempty(df_uniprot_res_date) #check that we still have row
            is_apo = is_apo = any(ismissing, df_uniprot_res_date.LIGANDS) && any(!ismissing, df_uniprot_res_date.LIGANDS)
            if is_apo #Check that we still have apo and holo conformation 
                return df_uniprot_res_date
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
For each uniprot accession we are adding information about the number of mutation and the missing residues, we calculate the RMSD 
between every pdb form the pdb accesion to compare there stucture and than cluster them in cluster base on there similarity

Take in input a dataframe for 1 uniprot and the path the the sift folder and the pdb folder 
Output a dataframe with all the added colonnes
"""
#function to be use in parallelization 
@everywhere function process_uniprot_df(df_uniprot, sift_folder, pdb_folder)
    #Filter the pdb 
    #df_uniprot_filter=filter_pdb(df_uniprot)
    #if df_uniprot_filter === nothing 
        #return missing
    #end
    df_uniprot_filter = count_mutation_pdb_uni(df_uniprot, sift_folder)

    rmsd_list, files, keep_indices = calculate_RMSD_uniprot(df_uniprot_filter, sift_folder, pdb_folder)

    if length(keep_indices) > 0
        df_cluster = cah(rmsd_list, files, keep_indices)
        select!(df_uniprot_filter, Not(["RES_END", "RES_BEG"]))
        df_result = leftjoin(df_uniprot_filter, df_cluster, on=[:PDB, :CHAIN])
        return df_result
    else
        return missing
    end
end

# Get all the information of each pdb 
"""
For each uniprot accession we are adding information about the number of mutation and the missing residues, we calculate the RMSD 
between every pdb form the pdb accesion to compare there stucture and than cluster them in cluster base on there similarity

Take in input the DataFrame from get_ligand_information
Output the DataFrame with new columns MUTATION MISSING_RESIDUES Cluster_1.5 Cluster_2.0 Cluster_4.0
"""
function clustering_each_uniprot_acc(df_final::DataFrame,sift_folder::String,pdb_folder::String)
    grouped_uniprots = DataFrame[g for g in groupby(df_final, :UNIPROT)]

    # Try only for 10
    grouped_uniprots_limited = grouped_uniprots[1:min(50, length(grouped_uniprots))]

    df_results = pmap(df -> process_uniprot_df(df, sift_folder, pdb_folder), grouped_uniprots_limited)
    df_completed = vcat(skipmissing(df_results)...)
    return df_completed
end

#Compare the Clustering 
"""
Function to analyse the different cutoff from the cluster and chose the more releatable one 

Take in input the DataFrame from clustering_each_uniprot_acc
Output a DataFrame to know for each Uniprot accession if we have the apo and holo form and for with cluster there where different
"""
function check_apo_holo_cluster(df_completed::DataFrame)
    check_cluster=DataFrame(UNIPROT=String[],HOLO_APO=Bool[],BEST_CUTOFF=Union{Float64,Missing}[]) # Create an empty DF

    for row in eachrow(df_completed)
        uniprot=row.UNIPROT # Check for every uniprot accession 
        #Check that we haven't already done this
        is_present= uniprot in check_cluster.UNIPROT
        if !is_present
            df_uniprot = select(filter(row -> row.UNIPROT == uniprot, df_completed), ["PDB", "LIGANDS", "Cluster_0.5",  "Cluster_0.75",  "Cluster_1.0",  "Cluster_1.25",  "Cluster_1.5","Cluster_2.0"])# To focus on the uniprot 
            
            #Look if we have both conformation for this uniprot
            is_apo = any(ismissing, df_uniprot.LIGANDS)
            if is_apo 
                # Filter files with and without ligands
                with_ligands = filter(row -> !ismissing(row.LIGANDS), df_uniprot)
                without_ligands = filter(row -> ismissing(row.LIGANDS), df_uniprot)
                result = 9
                # Compare clusters for each cutoff
                for cutoff in [0.5, 0.75, 1, 1.25, 1.5,2]
                    # Retrieve clusters for each cutoff
                    cluster_col = Symbol("Cluster_$(cutoff)")

                    # Compare files without ligands and with ligands
                    for file_with_ligands in eachrow(with_ligands)
                        for file_without_ligands in eachrow(without_ligands)
                            if !ismissing(file_with_ligands[cluster_col]) && 
                                !ismissing(file_without_ligands[cluster_col]) && 
                                file_with_ligands[cluster_col] != file_without_ligands[cluster_col]
                                # If the files are in different clusters, take the smaller value
                                if  cutoff<result
                                    result = cutoff
                                end
                            end
                        end
                    end
                end
                push!(check_cluster,(uniprot,true,result))
                
            else 
                push!(check_cluster,(uniprot,false,9))
            end 
        end
    end
    return check_cluster #Output the result 
end


#download sift Pfam and cath databases 
""" 
get the three csv that are already download and output them in a dataframe
"""
function download_file()
    ## Get Uniprot 
    sift_uniprot_mapping=get_uniprot_mapping()
    ## Get pfam
    sifts_file_path_pfam="pdb_chain_pfam.csv"
    sift_pfam_mapping=DataFrames.DataFrame(CSV.File(sifts_file_path_pfam,
            comment="#", missingstring=["", "None"])) #Output DF with PDB CHAIN SP_PRIMARY PFAM_ID COVERAGE

    ##Get CATH 
    sift_file_path_cath="pdb_chain_cath_uniprot.csv"
    sift_cath_mapping=DataFrames.DataFrame(CSV.File(sift_file_path_cath,
            comment="#", missingstring=["", "None"])) #Output DF with PDB CHAIN SP_PRIMARY CATH_ID

    return sift_uniprot_mapping, sift_pfam_mapping, sift_cath_mapping
end

#Take information of biolip file 
"""
take in input boolean to know if we already have the smaler file and if the save the result 
Output a dataframe with all the information we need from the BioLip databases 
"""
function biolip_preparation(create_file::Bool ,save::Bool)
    if create_file 
        ## Get biolip
        file_path_biolip="BioLiP.csv"
        df_biolip=CSV.File(file_path_biolip, delim='\t') |> DataFrame

        # Garder uniquement les 5 premières colonnes
        df_biolip_5_first_columns = select(df_biolip, 1:5)
        rename!(df_biolip_5_first_columns, :Column2=> "CHAIN")
        rename!(df_biolip_5_first_columns, :Column3=> "RESOLUTION")
        rename!(df_biolip_5_first_columns, :Column4=> "SITE")
        rename!(df_biolip_5_first_columns, :Column5=> "LIGAND")

        if save 
            # Sauvegarder en CSV
            CSV.write("Biolip_5_first_columns.csv", df_biolip_5_first_columns)
            println(first(df_biolip_5_first_columns,20))
            println(size(df_biolip_5_first_columns))
        end
    else 
        if save 
            ## Get 5 first columns of biolip
            file_path_biolip="Biolip_5_first_columns.csv"
            df_biolip_5_first_columns=DataFrames.DataFrame(CSV.File(file_path_biolip,
            comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
            println(first(df_biolip_5_first_columns,20))
            println(size(df_biolip_5_first_columns))
            rename!(df_biolip_5_first_columns, Symbol.(strip.(string.(names(df_biolip_5_first_columns))))) #Enlever les espaces dans le nom des colonnes 
        end 
    end 
    return df_biolip_5_first_columns
end


################  Function main ##################
""" 
main fonction to get all the information about all the pdb 
use PFAM CATH BioLiP SIFT 
and will after compare the pdb for each uniprot in clustering_each_uniprot_acc
Take in input Boolean to now f=if we do the first part and if we take a file save, take also the path to the folder 
where all the sift are save and all the pdb 
At the end we will have a csv created with all the information : pdb_information_details_final.csv
"""
function main(part_1::Bool,save::Bool,sift_folder::String,pdb_folder::String)
    sift_uniprot_mapping, sift_pfam_mapping, sift_cath_mapping = download_file()
    if part_1
        df_biolip_5_first_columns=biolip_preparation(false,false)

        ## Join PFAM and CATH
        sift_join_file=join_information(sift_pfam_mapping,sift_cath_mapping,sift_uniprot_mapping) #Output DF with PDB CHAIN RES_BEG RES_END SP_PRIMARY PFAM_NB CATH_NB

        ## Recuperer information pdb 
        df_pdb_reso_resi=get_pdb_information(sift_join_file) # Output DF with PDB CHAIN RES_BEG RES_END SP_PRIMARY PFAM_NB CATH_NB

        ##Recuperer information sur le ligand 
        df_final=get_ligand_information(df_biolip_5_first_columns,df_pdb_reso_resi)

        if save 
            CSV.write("pdb_information_details.csv", df_final)
        end
    end
    if save 
        # get back the final df 
        file_path_df_final="pdb_information_details_filter.csv"
        df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
        comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
    end 

    #Clustering
    df_completed=clustering_each_uniprot_acc(df_final,sift_folder,pdb_folder)
    println(first(df_completed,20))
    println(size(df_completed))

    
    # Save the result
    CSV.write("pdb_information_details_final_1.csv", df_completed)
    
    # Compare the cluster
    check_cluster=check_apo_holo_cluster(df_completed)
    println(first(check_cluster,20))
    value_counts = countmap(check_cluster.BEST_CUTOFF)  # Count the occurrences
    most_frequent_value = argmax(value_counts) 
    println("The cutoff that best separates the apo and holo form is : ", most_frequent_value)
    frequency = count(x -> !ismissing(x) && x == most_frequent_value, check_cluster.BEST_CUTOFF)
    println("She appears ", frequency, " time.")
    println(size(check_cluster))
    println("It separates properly : ",(frequency/size(check_cluster)[1])*100)
    CSV.write("pdb_apo_holo_1000.csv", check_cluster)
    
end

###########################################################################################################################################
########################################################## MAIN ############################################################################
#############################################################################################################################################

@info "START ! "
@show "Début de debut ", Dates.format(now(), "HH:MM:SS")

const sift_folder= abspath("/alpha/database/sift", "sift_files")
const pdb_folder= abspath("/alpha/database/pdb", "pdb_files")

main(false,true,sift_folder,pdb_folder)

@show "Début de fin", Dates.format(now(), "HH:MM:SS")

@info "END !"
