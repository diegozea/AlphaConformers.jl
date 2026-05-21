#!/home/julie.daniel/.julia/juliaup/julia-1.11.7+0.x64.linux.gnu/bin/julia

#=
#SBATCH --nodelist=node48
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --output=test_AlphaConformers.jl.o%j.out
#SBATCH --job-name=foldseek_test
=#
import Pkg

Pkg.activate("/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/Clean_AlphaConformers/scripts/update_MIToS_321")

Pkg.test("AlphaConformers")
