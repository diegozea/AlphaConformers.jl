#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=40
#SBATCH --output=create_databases-%j.out
=#
import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")
using Distributed
addprocs(12)  # Ajoute 4 processus parallèles

@everywhere using DataFrames, CSV
@everywhere using Dates
@everywhere using MIToS.PDB
@everywhere import MIToS
@everywhere using AlphaConformers
@everywhere using BioStructures
@everywhere using DataStructures
@everywhere using Statistics

@everywhere function check_rmsd_apo_holo(df_uniprot_res_date,pdb_folder::String)
    #Check that the apo and holo conformation have at least 1A of differences 
    holo = filter(row -> !ismissing(row.LIGANDS), df_uniprot_res_date)
    apo = filter(row -> ismissing(row.LIGANDS), df_uniprot_res_date)
    for row_holo in eachrow(holo)
        pdb_holo=row_holo.PDB
        pdbfilepath_holo=joinpath(pdb_folder,uppercase(pdb_holo)*".pdb" ) # Use the temporary file
        if isfile(pdbfilepath_holo) 
            #Get the right chain only
            residues_1ivo_holo = read(pdbfilepath_holo, PDBFile)
            chain_residues_pdb_holo = MIToS.PDB.residuesdict(residues_1ivo_holo; model="1", chain=String(row_holo.CHAIN),  group="ATOM") #Returns dictionary with position => residues details 
            rmsd_list = []
            for row_apo in eachrow(apo)
                pdb_apo=row_apo.PDB
                pdbfilepath_apo=joinpath(pdb_folder,uppercase(pdb_apo)*".pdb" ) # Use the temporary file
                if isfile(pdbfilepath_apo) 
                    #Get the right chain only
                    residues_1ivo_apo = read(pdbfilepath_apo, PDBFile)
                    chain_residues_pdb_apo =MIToS.PDB.residuesdict(residues_1ivo_apo; model="1", chain=String(row_apo.CHAIN),  group="ATOM") #Returns dictionary with position => residues details 
                    #get the alignment 
                    common_keys = Set(keys(chain_residues_pdb_holo)) ∩ Set(keys(chain_residues_pdb_apo))
                    positions = SortedSet(parse.(Int,common_keys))  # Retrieves common positions (same key in all files)
                    # Extracting residues for each file as an ordered list
                    residues_holo = [chain_residues_pdb_holo[string(pos)] for pos in positions] 
                    residues_apo = [chain_residues_pdb_apo[string(pos)] for pos in positions]
                    rmsd_val=missing
                    try  
                        # Superposition and calculation of RMSD
                        _, _,rmsd_val = superimpose(residues_holo, residues_apo)
                    catch e 
                        @show common_keys
                        @warn "❌ Error while superimposing structures "
                        # we will not annalyse those pdb 
                    end 
                    push!(rmsd_list, rmsd_val)
                else 
                    @warn "The PDB file of $pdb_apo doesn't exist"
                    filter!(row -> !(row.PDB == pdb_apo), df_uniprot_res_date)

                end
            end
            if mean(rmsd_list)<1
                filter!(row -> !(row.PDB == pdb_holo && row.CHAIN == String(row_holo.CHAIN)), df_uniprot_res_date)

            end
        else 
            @warn "The PDB file of $pdb_holo doesn't exist"
            filter!(row -> !(row.PDB == pdb_holo), df_uniprot_res_date)
        end 
    end 
    is_apo = is_apo = any(ismissing, df_uniprot_res_date.LIGANDS) && any(!ismissing, df_uniprot_res_date.LIGANDS)
    if is_apo
        return df_uniprot_res_date
    else 
        return nothing
    end
end

@everywhere function filter_pdb(df_uniprot,pdb_folder::String)
    # Check that apo and holo are in different cluster 
    with_ligands = filter(row -> !ismissing(row.LIGANDS), df_uniprot)
    without_ligands = filter(row -> ismissing(row.LIGANDS), df_uniprot)
    result = 0
    for cutoff in [0.5, 0.75, 1.0, 1.25, 1.5,2]
        cluster_col = Symbol("Cluster_$(cutoff)")
        for file_with_ligands in eachrow(with_ligands)
            for file_without_ligands in eachrow(without_ligands)
                cluster_with = file_with_ligands[cluster_col]
                cluster_without = file_without_ligands[cluster_col]
                if !ismissing(cluster_with) && !ismissing(cluster_without) && cluster_with != cluster_without
                    # If the files are in different clusters, take the smaller value
                    if  cutoff>result
                        result = cutoff
                    end 
                end
            end
        end
    end
    if result != 0
        #check that they have more than 1A of differences 
        if result >=1.0
            cutoff =1.0
            cluster_col = Symbol("Cluster_$(cutoff)")
           # Extraire les clusters (non-missing) pour cutoff 1.0
            clusters_with = Set(row[cluster_col] for row in eachrow(with_ligands) if !ismissing(row[cluster_col]))
            clusters_without = Set(row[cluster_col] for row in eachrow(without_ligands) if !ismissing(row[cluster_col]))

            # Garder les lignes seulement si les clusters sont disjoints
            if isempty(intersect(clusters_with, clusters_without))
                df_final = vcat(with_ligands, without_ligands)
            else
                df_final = DataFrame()  # Vide si clusters en commun
            end
            @show df_final
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

function create_database(df_final::DataFrame,pdb_folder:: String)
    grouped_uniprots = groupby(df_final, :UNIPROT)
    chunks = Iterators.partition(grouped_uniprots, 1000)  # 1000 groupes à la fois

    results = pmap(chunk -> begin
        tmp = Vector{DataFrame}()
        for group in chunk
            res = filter_pdb(group,pdb_folder)
            if res !== nothing
                push!(tmp, res)
            end
        end
        tmp
    end, collect(chunks))

    return vcat(Iterators.flatten(results)...)
end
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
        return date_info
    else
        return nothing
    end
    
end

function foldseek_similar_pdb(df_final::DataFrame,FOLDSEEK_DB::String,pdb_folder::String)
    database=DataFrame(UNIPROT=String[],PDB=String[],CHAIN=String[],NB_SIMILAR_PROT=Int64[]) # Create an empty DF
    df_merged=DataFrame()
    for row in eachrow(df_final)
        uniprot=row.UNIPROT
        pdb_file=row.PDB
        chain=row.CHAIN
        println(pdb_file)
        pdbfilepath=joinpath(pdb_folder,uppercase(pdb_file)*".pdb" ) # Use the temporary file
        if isfile(pdbfilepath) #Check that the file is not already downloaded
            @info "Run Foldseek"
            out_file=AlphaConformers.foldseek_search(pdbfilepath,db_path=FOLDSEEK_DB)
            df_result=AlphaConformers.read_foldseek_search_results(out_file)
            #Filter the evalue to keep the one < 1e-3
            println(first(df_result,5))
            df_result_evalue = filter(row -> row.evalue < 1e-5, df_result)
            count=0
            for row in eachrow(df_result_evalue)
                #Check that the publication date is before AF2 training
                target=String(split(row.target,".")[1])
                code_pdb, date_info, method_info, res_info = get_header_pdb(target)
                if Date(date_info) < Date("2018-04-30")
                    count +=1
                end
            end
            push!(database,(uniprot,pdb_file,chain,count))
        else 
            @warn "The PDB file of $pdb_file doesn't exist"
        end
        
    end
    println(database)
    df_merged = innerjoin(database, df_final, on=[:UNIPROT, :PDB, :CHAIN])
    return df_merged
end

function use_usalign(df_result_evalue,pdbfilepath)
    for (i, row2) in enumerate(eachrow(df_result_evalue))
        if i != 1
            row2.target
            aln_data = mktempdir() do tmp_folder
                aln_data = nothing
                try 
                    pdb, chain = split(row2.target,"_")
                    @info "Downloading $pdb"
                    pdb_file = joinpath(tmp_folder, "$pdb.gz")
                    MIToS.PDB.downloadpdb(String(split(pdb,".")[1]); format=MIToS.PDB.PDBFile, filename=pdb_file)
                    aln_data = usalign(pdbfilepath, pdb_file)
                    
                catch e 
                    @warn "Didn't succeed to download the pdb"
                    aln_date=nothing
                end
                return aln_data
            end
            if aln_data === nothing
                return false
            end
            @show first(aln_data,5)
            row=first(aln_data)
            if row.TM1 > 0.5 || row.TM2 > 0.5
                return  false
            end
        end
    end
    return true
end
@info "Start"
const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"
pdb_folder= abspath("/alpha/database/pdb", "pdb_files")

#Get the file with the result 
file_path_df_final="pdb_information_details_final_mutation_cluster_reformatted.csv"
df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(df_final,20))
println(size(df_final))

filtering = true

if filtering 
    df_final=create_database(df_final,pdb_folder)
    CSV.write("pdb_information_details_final_mutation_cluster_reformatted_filter.csv", filter(!isnothing, df_final)) 
end
println(first(df_final,20))
println(size(df_final))


if df_final !== nothing 
    df_merged=foldseek_similar_pdb(df_final,FOLDSEEK_DB,pdb_folder)
    println(first(df_merged,20))
    println(size(df_merged))
    CSV.write("pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek.csv",df_merged)
end

if df_merged !== nothing
    database = filter(row -> row.NB_SIMILAR_PROT < 50, df_merged)
    @show size(database)
    df=DataFrame(UNIPROT=String[],PDB=String[],CHAIN=String[]) # Create an empty DF
    df_last=DataFrame()
    for row in eachrow(database)
        uniprot=row.UNIPROT
        pdb_file=row.PDB
        chain=row.CHAIN
        println(pdb_file)
        pdbfilepath=joinpath(pdb_folder,uppercase(pdb_file)*".pdb" ) # Use the temporary file
        if isfile(pdbfilepath) #Check that the file is not already downloaded
            @info "Run Foldseek"
            out_file=AlphaConformers.foldseek_search(pdbfilepath,db_path=FOLDSEEK_DB)
            df_result=AlphaConformers.read_foldseek_search_results(out_file)
            #Filter the evalue to keep the one < 1e-3
            println(first(df_result,5))
            df_result_evalue = filter(row -> row.evalue < 1e-5, df_result)
            filtered_df = filter(row -> begin
                target = split(row.target, ".")[1]
                code_pdb, date_info, _, _ = get_header_pdb(target)
                row.evalue < 1e-5 && Date(date_info) < Date("2018-04-30")
            end, df_result_evalue)
            check=use_usalign(filtered_df,pdbfilepath)
            if check
                push!(df,(uniprot,pdb_file,chain))
            end
        else 
            @warn "The PDB file of $pdb_file doesn't exist"
        end
    end
    df_last = innerjoin(database, df, on=[:UNIPROT, :PDB, :CHAIN])
    CSV.write("pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek_final.csv", df_last) 
end


filtering = true

if filtering 
    df_final=create_database(df_final,pdb_folder)
    println(first(df_final,20))
    println(size(df_final))
    CSV.write("pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek_final_1.csv", filter(!isnothing, df_final)) 
end


