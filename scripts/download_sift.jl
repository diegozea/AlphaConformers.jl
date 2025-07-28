#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#SBATCH --nodelist=node48
#SBATCH --time=300:00:00
#SBATCH --output=download_sift.jl.%j.out
=#
import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")

# Load necessary packages on all workers
using BioStructures
using MIToS.PDB
using MIToS.SIFTS
using DataFrames
import CSV

####################################################### MAIN ########################################################
"""
This script downloads SIFT XML files from the EBI FTP server and saves them in a specified directory.
"""

@info "START"
# Define the directory where SIFT files will be saved
sift_directory = abspath("/alpha/database", "sift")
isdir(sift_directory) || mkdir(sift_directory)
output_directory = abspath(sift_directory, "sift_files")
isdir(output_directory) || mkdir(output_directory)

try
    run(`wget -r -np -nH --cut-dirs=5 -P $output_directory ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/`)
    @info "✅ Tous les fichiers XML SIFTS ont été téléchargés dans $output_directory"
catch e
    @warn "❌ Erreur pendant le téléchargement via wget: $e"
end

@info "END !"
####################################################### END ########################################################