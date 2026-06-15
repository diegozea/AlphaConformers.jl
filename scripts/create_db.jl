#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l ncpus=20
#PBS -l mem=100g
#PBS -l host=node48
#PBS -l walltime=300:00:00
#PBS -j oe
=#

# Set the working directory
cd("/alpha/database/pdb/")

const foldseek = "/store/EQUIPES/AMIG/MEMBERS/carla.martins/foldseek/foldseek/bin/foldseek"

# Run the foldseek createdb command
run(`$foldseek createdb pdb_files fullpdb --threads 16`)

# Run the foldseek createindex command
run(`$foldseek createindex fullpdb tmp --threads 16`)

println("THE END!")
