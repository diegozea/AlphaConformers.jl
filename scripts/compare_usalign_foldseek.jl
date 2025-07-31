#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#SBATCH --nodelist=node48
#SBATCH --time=10:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --output=compare_usalign_foldseek.jl.o%j.out
=#

import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/") 

# Load necessary packages
using CSV
using DataFrames
using StatsPlots
using Plots
using Measures
using Statistics
using AlphaConformers
using MIToS
using Foldseek_jll

########################################## MAIN #######################################################"
"""
This code compare the tm-score of usalign and foldseek for the same pdbs

Input : 
- pdb1 : path to the first pdb 
- pdb2 : path to the second pdb
Output : 
- Compare the two TM-score output by Foldseek and USalign 

"""
########################## Information to fill #################################

# Path to the 2pdb to compare - need to be save in the directory
pdb1 = "1HW1_B_1.pdb"
pdb2 = "1H9G_A_1.pdb"

################################################################################

# --- USAlign ---
usalign_df=usalign(pdb1, pdb2)
@show usalign_df

tm_usalign = usalign_df[!, :TM1][1]  

println("USAlign TM-score: ", tm_usalign)

# --- Foldseek ---
run(`$(Foldseek_jll.foldseek()) easy-search $pdb1 $FOLDSEEK_DB foldseek_result.m8 /store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/tmp --format-output query,target,alntmscore`)

foldseek_df = CSV.read("foldseek_result.m8", DataFrame, delim='\t', header=[
    "query", "target", "tmscore"
])

tm_foldseek = foldseek_df[!, :tmscore][1]

println("Foldseek TM-score: ", tm_foldseek)

# Comparaison
if abs(tm_usalign - tm_foldseek) < 1e-4
    println("Les TM-score sont identiques (ou très proches) !")
else
    println("Les TM-score sont différents.")
end
###################################################### END ################################################################