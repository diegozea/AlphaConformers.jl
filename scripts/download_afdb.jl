#!/stockage/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l ncpus=40
#PBS -l mem=200g
#PBS -l host=node48
#PBS -l walltime=999:99:99
#PBS -j oe
=#

cd("/alpha/database/afdb")

run(`/stockage/EQUIPES/AMIG/MEMBERS/carla.martins/foldseek/foldseek/bin/foldseek databases --compressed 1 -v 3 --threads 40 Alphafold/UniProt afdb_up tmp`)
