#!/bin/bash
#SBATCH --nodelist=node48
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --output=test_AlphaConformers.jl.o%j.out
#SBATCH --job-name=foldseek_test
#SBATCH --chdir=/store/EQUIPES/AMIG/MEMBERS/julie.daniel/Clean_AlphaConformers

#=
/home/julie.daniel/.julia/juliaup/julia-1.11.7+0.x64.linux.gnu/bin/julia \
  --project=/store/EQUIPES/AMIG/MEMBERS/julie.daniel/Clean_AlphaConformers \
  -e 'using TestItems, TestItemRunner; include("test/foldseek.jl"); @run_package_tests verbose=true'
=#