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
Pkg.status("MIToS")
using MIToS
using MIToS.PDB
using DataFrames
import CSV
using Revise
using AlphaConformers


const PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data"
cd(PATH)

const PDB_FOLDER = nothing

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"

# download the dataset 
file_path_df_final="/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/selected_examples.csv"
df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(df_final,20))
println(size(df_final))
"""
info_pdb = DataFrame(
    PDB_apo = String[], 
    CHAIN_apo = String[], 
    INDEX_apo = Union{Missing, String}[],  
    PDB_holo = String[], 
    CHAIN_holo = String[], 
    INDEX_holo = Union{Missing, String}[] 
)
for row in eachrow(df_final)
    apo=row.apo_id
    apo_info=split(apo,r"[_-]")
    if length(apo_info)>2
        apo_pdb=apo_info[1]
        apo_chain=apo_info[3]
        apo_index=apo_info[2]
    else
        apo_pdb=apo_info[1]
        apo_chain=apo_info[2]
        apo_index=missing
    end
    holo=row.holo_id
    holo_info=split(holo,r"[_-]")
    if length(holo_info)>2
        holo_pdb=holo_info[1]
        holo_chain=holo_info[3]
        holo_index=holo_info[2]
    else
        holo_pdb=holo_info[1]
        holo_chain=holo_info[2]
        holo_index=missing
        
    end
    push!(info_pdb,(apo_pdb,apo_chain,apo_index,holo_pdb,holo_chain,holo_index))
    println(holo_pdb)
    apo_path=joinpath(PATH,apo_pdb*".pdb.gz" )
    MIToS.PDB.downloadpdb(String(apo_pdb), format = PDBFile,filename=apo_path) # Sinon ne reconnais pas le package --> weird   
    
    if ismissing(apo_index)
        apo_index="1"
    else 
        apo_index=String(apo_index)
    end
    
    chain = AlphaConformers._read_pdb_chain(apo_path, String(apo_chain))
    if chain !== nothing
        apo_chain_path=joinpath(PATH,apo_pdb*"_"*apo_chain*".pdb.gz" )
        MIToS.PDB.write(apo_chain_path, chain, MIToS.PDB.PDBFile)
    end
    holo_path=joinpath(PATH,holo_pdb*".pdb.gz" )
    MIToS.PDB.downloadpdb(String(holo_pdb), format = PDBFile,filename=holo_path) # Sinon ne reconnais pas le package --> weird   
    
    if ismissing(holo_index)
        holo_index="1"
    else 
        holo_index=String(holo_index)
    end
    

    chain = AlphaConformers._read_pdb_chain(holo_path, String(holo_chain))
    if chain !== nothing
        holo_chain_path=joinpath(PATH,holo_pdb*"_"*holo_chain*".pdb.gz" )
        MIToS.PDB.write(holo_chain_path, chain, MIToS.PDB.PDBFile)
    end
end

println(info_pdb)
CSV.write("info_dev_set.csv", info_pdb)
"""
file_path="info_dev_set.csv"
info_pdb=DataFrames.DataFrame(CSV.File(file_path,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(info_pdb,20))
println(size(info_pdb))

#Run AlphaConformers with apo form in template 
global index=0
for row in eachrow(info_pdb)
    if index == 1
        apo_pdb=row.PDB_apo
        apo_chain=row.CHAIN_apo
        println(apo_pdb)
        REF_PDB = joinpath(PATH, apo_pdb*"_"*apo_chain*".pdb.gz")
        println(REF_PDB)
        output_dir = joinpath(PATH, apo_pdb)
        if isdir(output_dir)
            rm(output_dir; recursive=true, force=true)
        end
        mkdir(output_dir)
        println(output_dir)
        AlphaConformers.alphaconformers(REF_PDB, FOLDSEEK_DB, ALPHAFOLD_DB, PDB_FOLDER, output_dir)
        AlphaConformers.run_alphafold(output_dir, colabfold_path=COLABFOLD_PATH)
    end
    global index=index+1
end
