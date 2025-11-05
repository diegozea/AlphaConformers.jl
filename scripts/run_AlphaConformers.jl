#!/home/julie.daniel/.julia/juliaup/julia-1.11.7+0.x64.linux.gnu/bin/julia


#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --output=run_AlphaConformers.jl.o%j.out
=#


import Pkg


Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/Clean_AlphaConformers/scripts/update_MIToS_321")
Pkg.status("MIToS")

using AlphaConformers
@show pathof(AlphaConformers)
#Pkg.develop(path="/store/EQUIPES/AMIG/MEMBERS/julie.daniel/Clean_AlphaConformers")




# Load necessary packages 
using MIToS
using MIToS.PDB
using DataFrames
import CSV
using Revise
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
const PDB_FOLDER = "/alpha/database/pdb/mmcif_files"
# Path to the Foldseek database
const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb_mmcif_files"
# Path to the AlphaFold database
const ALPHAFOLD_DB = "/alpha/database/afdb_v6/fullafdb_v6"
#const ALPHAFOLD_DB = nothing
# Path to the ColabFold script
const SIF_PATH=expanduser("/store/EQUIPES/AMIG/SCRIPTS/sif_images/ColabFold_AF2_1-5-5/colabfold-1.5.5-cuda12.2.2.sif")
const CACHE_DIR="/store/EQUIPES/AMIG/SCRIPTS/sif_images/ColabFold_AF2_1-5-5/cache"
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
output_dir = joinpath(PATH, apo_pdb*"_AlphaConformer_new_bdd")


if isdir(output_dir)
    rm(output_dir; recursive=true, force=true)
end
mkdir(output_dir)
println(output_dir)

#Run AlphaConformers

AlphaConformers.alphaconformers(REF_PDB, PDB_FOLDER, output_dir; db=db,evalue_cutoff=NaN, cutoff=1.0)

AlphaConformers.run_alphafold_one_run(output_dir, SIF_PATH, CACHE_DIR)   
         

@show "End"
########################################### End ###########################################################