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

function _read_pdb_chain_model(file::String, chain_code::String,model_code::String)
    # Read the chain specified by chain_code from the PDB file
    # occupancyfilter is needed to avoid the duplicated residue warnings with TMalign
    try
        # read the whole file
        res = read(file, MIToS.PDB.PDBFile, onlyheavy=true, occupancyfilter=true)
        # note that auth chains can be lower case, so test the lowercase one if 
        # the uppercase one is not found. For example, 7ADD has lowercase chains.
        chains = Set{String}(r.id.chain for r in res)
        model = Set{String}(r.id.model for r in res) #If no model specify need to put chain_code="1"
        if chain_code in chains && model_code in model
            MIToS.PDB.residues(res, model_code, chain_code)
            
        else
            lowercase_chain_code = lowercase(chain_code)
            if lowercase_chain_code in chains && model_code in model
                return MIToS.PDB.residues(res, model_code, lowercase_chain_code)
            end
            @error "The chain $chain_code or the model $model_code was not found in the PDB file $file"
            nothing
        end
    catch err
        @error "Error reading the PDB file $file: $err"
        nothing
    end
end

const PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data"
cd(PATH)

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"
#const ALPHAFOLD_DB = nothing

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
        apo_index="1"
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
        holo_index="1"
        
    end
    push!(info_pdb,(apo_pdb,apo_chain,apo_index,holo_pdb,holo_chain,holo_index))
    println(holo_pdb)

    apo_path=joinpath(PATH,apo_pdb*".pdb.gz" )
    MIToS.PDB.downloadpdb(String(apo_pdb), format = PDBFile,filename=apo_path) # Sinon ne reconnais pas le package --> weird   
    
    chain = _read_pdb_chain_model(apo_path, String(apo_chain),String(apo_index))
    if chain !== nothing
        apo_chain_path=joinpath(PATH,apo_pdb*"_"*apo_chain*"_"*apo_index*".pdb.gz" )
        MIToS.PDB.write(apo_chain_path, chain, MIToS.PDB.PDBFile)
    end
    
    holo_path=joinpath(PATH,holo_pdb*".pdb.gz" )
    MIToS.PDB.downloadpdb(String(holo_pdb), format = PDBFile,filename=holo_path) # Sinon ne reconnais pas le package --> weird   

    chain = _read_pdb_chain_model(holo_path, String(holo_chain),String(holo_index))
    if chain !== nothing
        holo_chain_path=joinpath(PATH,holo_pdb*"_"*holo_chain*"_"*holo_index*".pdb.gz" )
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
    if index == 4
        apo_pdb=row.PDB_apo
        apo_chain=row.CHAIN_apo
        apo_model=row.INDEX_apo
        println(typeof(apo_model))
        println(apo_pdb)
        REF_PDB = joinpath(PATH, apo_pdb*"_"*apo_chain*"_"*string(apo_model)*".pdb.gz")
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