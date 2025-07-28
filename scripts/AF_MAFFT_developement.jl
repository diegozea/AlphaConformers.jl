#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10
#SBATCH --output=AF_MAFFT_developement.jl%j.out
=#
import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")

# Load necessary packages
using FilePathsBase, Glob
using AlphaConformers
import MAFFT_jll
using MIToS.MSA

######################################################## MAIN  #####################################################
"""
This scripts change the MSA output by Foldseek by using MAFFT
Take in input the ouput of AlphaConformer extract everything 
And change only the MSA 
With MAFFT we align by sequence and not structure 

Input : 
- PATH : path to the folder with all the output 
- COLABFOLD_PATH : PATH to the colab fold to run  AF2
- source_root : PATH to the output of AlphaConformer
- dest_root : PATH to the new folder created to apply MAFFT
Output : 
Same file organization than ALphaConformer with MSA align with MAFFT 
HAve the AF2 output 

This code was use as a test 
Will need to change how it work to be more efficent
"""

########################## Information to fill #################################
#Path where all the result are
const PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data"
#Path to colabfold
const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
# Path to get the MSA and template
source_root = joinpath(PATH,"1AKZ_No_AFDB")
#Path to save the result
dest_root = joinpath(PATH,"1AKZ_No_AFDB_MAFFT")
###################################################################################

cd(PATH)

# create the output folder
mkpath(dest_root)

# For all the cluster 
for cluster_dir in glob("cluster_*", source_root)
    src_cluster_path = joinpath(source_root, cluster_dir)
    dest_cluster_path = joinpath(dest_root, basename(cluster_dir))
    @show cluster_dir
    @show src_cluster_path
    @show dest_cluster_path
    # Create the output folder with the cluster
    mkpath(dest_cluster_path)

    # 1. Copie the template folder
    src_templates = joinpath(src_cluster_path, "templates")
    dest_templates = joinpath(dest_cluster_path, "templates")
    mkpath(dest_templates)
    for pdb_file in readdir(src_templates)
        cp(joinpath(src_templates, pdb_file), joinpath(dest_templates, pdb_file); force=true)
    end
    # 2.Apply MAFFT on the MSA
    src_seq_file = joinpath(src_cluster_path, "sequences.a3m")
    dest_seq_file = joinpath(dest_cluster_path, "sequences.a3m")

    if isfile(src_seq_file)
        run(pipeline(`$(MAFFT_jll.mafft()) --quiet $src_seq_file`, stdout=dest_seq_file))
    else
        @warn "Fichier sequences.a3m manquant dans $src_cluster_path"
    end
    #Take the output to save in the write shape 
    msa = read(dest_seq_file, A3M)
    output_file=adjustreference!(msa)
    write(dest_seq_file,output_file,A3M)
end

#Run AlphaFold
AlphaConformers.run_alphafold(dest_root, colabfold_path=COLABFOLD_PATH)

############################################################ END ##############################################################
