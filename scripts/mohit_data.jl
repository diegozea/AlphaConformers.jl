#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#PBS -l host=node48
#PBS -l walltime=900:00:00
#PBS -l mem=100gb
#PBS -l ncpus=40
#PBS -j oe
=#

using Revise
using CSV, DataFrames
using AlphaConformers

const PATH = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/mohit"
cd(PATH)

const PDB_FOLDER = nothing

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"

output_dir = joinpath(PATH, "first_run")
if isdir(output_dir)
    rm(output_dir; recursive=true, force=true)
end
mkdir(output_dir)

const REF_PDB = joinpath(PATH, "P51993_35_359.pdb")

AlphaConformers.alphaconformers(REF_PDB, FOLDSEEK_DB, ALPHAFOLD_DB, PDB_FOLDER, output_dir)
AlphaConformers.run_alphafold(output_dir, colabfold_path=COLABFOLD_PATH)
