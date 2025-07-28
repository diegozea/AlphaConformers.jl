#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --output=AlphaConformers_development_test.jl.o%j.out
=#


import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")
Pkg.status("MIToS")

# Load necessary packages
using MIToS
using MIToS.PDB
using DataFrames
import CSV
using Revise
using AlphaConformers

################################################## Functions ######################################################
"""
Extract the chain from a PDB file and save it to a new file.
This function reads a PDB file, extracts the specified chain, and writes it to a new PDB file.

"""
# Function to read a PDB file and extract a specific chain
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

# Function to create a development set
function create_dev_set()
    # download the original dataset 
    file_path_df_final="/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/selected_examples.csv"
    df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
    comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
    println(first(df_final,20))
    println(size(df_final))
    
    #Create the DF
    info_pdb = DataFrame(
        PDB_apo = String[], 
        CHAIN_apo = String[], 
        INDEX_apo = Union{Missing, String}[],  
        PDB_holo = String[], 
        CHAIN_holo = String[], 
        INDEX_holo = Union{Missing, String}[] 
    )

    #For each pdb we want 
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
        #save the information
        push!(info_pdb,(apo_pdb,apo_chain,apo_index,holo_pdb,holo_chain,holo_index))
        println(holo_pdb)

        #Download the apo and holo pdb files with the chain and index
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

    #Save the information in a csv file
    CSV.write("info_dev_set.csv", info_pdb)

end

################################################## MAIN ######################################################
"""
Function to executes AlphaConformer on a set of PDB files

Input :
-Path: Path to the main directory containing the apo and holo files
-PDB_FOLDER: Path to the PDB files
-FOLDSEEK_DB: Path to the Foldseek database
-ALPHAFOLD_DB: Path to the AlphaFold database
-COLABFOLD_PATH: Path to the ColabFold script
-db: BDD use for Foldseek
Output:
- Folder for each PDB file containing the results of AlphaConformer
- Function returns cluster and each have the result of AlphaFold

Need to wait for AlphaFold result 
Take around 40min for each pdb for AlphaConformer then wait for AlphaFold to execute each cluster 
Time of total execution depend on the number of cluster 
Can change the parameter of ALphaConformer to reduce the number of cluster
    - db, evalue_cutoff, cutoff
"""
########################## Information to fill #################################
# Path to the main directory containing the apo and holo files and where to save the result
const PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/"
cd(PATH)
# Path to the PDB files
const PDB_FOLDER = "/alpha/database/pdb/pdb_files"
# Path to the Foldseek database
const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"
# Path to the AlphaFold database
const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"
#const ALPHAFOLD_DB = nothing
# Path to the ColabFold script
const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
#BDD use for Foldseek
db=[FOLDSEEK_DB,ALPHAFOLD_DB]
####################################################################################

#Create the development set
if !isfile("info_dev_set.csv")
    create_dev_set()
end
file_path="info_dev_set.csv"
info_pdb=DataFrames.DataFrame(CSV.File(file_path,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(info_pdb,20))
println(size(info_pdb))

#Run AlphaConformers with apo form in template 
global index=0
for row in eachrow(info_pdb)
    
    #Get the apo pdb file
    apo_pdb=row.PDB_apo
    apo_chain=row.CHAIN_apo
    apo_model=row.INDEX_apo #if model specified
    @show apo_pdb
    filename = string(apo_pdb, "_", apo_chain, "_", apo_model, ".pdb.gz")
    REF_PDB = joinpath(PATH, filename) #Get the query path
    @show REF_PDB

    #Create the output directory
    #output_dir = joinpath(PATH, apo_pdb*"_No_AFDB")
    #output_dir = joinpath(PATH, apo_pdb)
    output_dir = joinpath(PATH, "Devset_test/"*apo_pdb*"_AlphaConformer_Hobohm")
    if isdir(output_dir)
        rm(output_dir; recursive=true, force=true)
    end
    mkdir(output_dir)
    println(output_dir)
    
    #Run AlphaConformers
    try 
        AlphaConformers.alphaconformers(REF_PDB, PDB_FOLDER, output_dir; db=db, cutoff=1.0)
        AlphaConformers.run_alphafold(output_dir, colabfold_path=COLABFOLD_PATH)   
    catch e
        @warn "Error for $apo_pdb : ", e
    end         
    
    global index=index+1 #can be use if we want to run only for few pdbs
end
@show "End"
########################################### End ###########################################################