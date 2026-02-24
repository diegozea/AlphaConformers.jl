#!/stockage/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l host=node48
#PBS -l walltime=300:00:00
#PBS -j oe
=#

using BioStructures

output_directory = abspath("/alpha/database/pdb", "pdb_files")
isdir(output_directory) || mkdir(output_directory)

downloadentirepdb(dir=output_directory)

# delete all the PDB files that do not contain the string "ATOM"
# diego.zea@node48:/alpha/database/pdb/pdb_files$ grep -L "ATOM" * | xargs rm -f

println("THE END!")