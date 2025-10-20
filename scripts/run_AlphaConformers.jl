#!/home/julie.daniel/.julia/juliaup/julia-1.11.7+0.x64.linux.gnu/bin/julia


#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --output=run_AlphaConformers.jl.o%j.out
=#


import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/Clean_AlphaConformers/scripts/update")
Pkg.status("MIToS")

# Load necessary packages 
using MIToS
using MIToS.PDB
using DataFrames
import CSV
using Revise
using AlphaConformers
using Glob
using Statistics

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
const PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/Clean_AlphaConformers/data/"
cd(PATH)
# Path to the PDB files
const PDB_FOLDER = "/alpha/database/pdb/pdb_files"
# Path to the Foldseek database
const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"
# Path to the AlphaFold database
const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"
#const ALPHAFOLD_DB = nothing
# Path to the ColabFold script
const COLABFOLD_PATH = "/opt/alphafold/scripts/runcolabfold.py"
#BDD use for Foldseek
db=[FOLDSEEK_DB,ALPHAFOLD_DB]
####################################################################################

#Run AlphaConformers with apo form in template 

apo_pdb="6r17"
apo_chain="C"

filename = string(apo_pdb, "_", apo_chain, ".pdb")
REF_PDB = joinpath(PATH, filename) #Get the query path
@show REF_PDB

#Create the output directory
#output_dir = joinpath(PATH, apo_pdb*"_No_AFDB")
#output_dir = joinpath(PATH, apo_pdb)
output_dir = joinpath(PATH, apo_pdb*"_AlphaConformer")


if isdir(output_dir)
    rm(output_dir; recursive=true, force=true)
end
mkdir(output_dir)
println(output_dir)

#Run AlphaConformers

AlphaConformers.alphaconformers(REF_PDB, PDB_FOLDER, output_dir; db=db,evalue_cutoff=NaN, cutoff=1.0)

#=
## Code to have only one run AF2 ##
#Create a folder to gather all the a3m files
a3m_folder = joinpath(output_dir, "All_a3m")
if isdir(a3m_folder)
    rm(a3m_folder; recursive=true, force=true)
end
mkdir(a3m_folder)

#Create a folder to gather all the a3m files
template_all_folder = joinpath(output_dir, "All_templates")
if isdir(template_all_folder)
    rm(template_all_folder; recursive=true, force=true)
end
mkdir(template_all_folder)

#Get the cluster folders
cluster_folders = filter!(
        dir -> occursin("cluster_", dir), 
        readdir(clusters_folder, join=true))
@show cluster_folders
#Copy the a3m files in the new folder
for folder in cluster_folders
    a3m_file = joinpath(folder, "sequences.a3m")
    @show a3m_file
    if isfile()
        cp(a3m_file, joinpath(a3m_folder, "msa_"*folder[end]*".a3m"))
    else
        @warn "No a3m file found in $folder"
    end
    template_folder = joinpath(folder, "templates")
    @show template_folder
    if isfolder()
        for file in readdir(template_folder)
            if !isfile(joinpath(template_all_folder, file))
                cp(joinpath(template_folder, file), joinpath(template_all_folder, file))
            end
        end
    else
        @warn "No template folder found in $folder"
    end
end

AlphaConformers.run_alphafold_one_run(output_dir, colabfold_path=COLABFOLD_PATH)   
=#
###################################################

#AlphaConformers.run_alphafold(output_dir, colabfold_path=COLABFOLD_PATH)   
         

@show "End"
########################################### End ###########################################################