#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l host=node48
#PBS -l walltime=300:00:00
#PBS -j oe
=#
import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")

using BioStructures
using MIToS.PDB
using MIToS.SIFTS
using DataFrames
import CSV
@info "START"
sift_directory = abspath("/alpha/database", "sift")
isdir(sift_directory) || mkdir(sift_directory)
output_directory = abspath(sift_directory, "sift_files")
isdir(output_directory) || mkdir(output_directory)

output_directory_pdb = abspath("/alpha/database/pdb", "pdb_files")

#get the pdb we want 
file_path_df_final="pdb_information_details.csv"
df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND

try
    run(`wget -r -np -nH --cut-dirs=5 -P $output_directory ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/`)
    @info "✅ Tous les fichiers XML SIFTS ont été téléchargés dans $output_directory"
catch e
    @warn "❌ Erreur pendant le téléchargement via wget: $e"
end

@info "END !"