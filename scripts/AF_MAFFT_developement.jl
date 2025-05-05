#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#PBS -l host=node48
#PBS -l walltime=900:00:00
#PBS -l mem=100gb
#PBS -l ncpus=40
#PBS -j oe
=#

import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
Pkg.add("FilePathsBase")

using FilePathsBase, Glob
using AlphaConformers
import MAFFT_jll
using MIToS.MSA

const PATH = "/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data"
cd(PATH)
const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
# Chemins
source_root = joinpath(PATH,"1AKZ_No_AFDB")
dest_root = joinpath(PATH,"1AKZ_No_AFDB_MAFFT")
@show source_root
@show dest_root

# Création de la racine de sortie si elle n'existe pas
mkpath(dest_root)

# Parcours des clusters
for cluster_dir in glob("cluster_*", source_root)
    src_cluster_path = joinpath(source_root, cluster_dir)
    dest_cluster_path = joinpath(dest_root, basename(cluster_dir))
    @show cluster_dir
    @show src_cluster_path
    @show dest_cluster_path
    # Crée le dossier de destination pour le cluster
    mkpath(dest_cluster_path)

    # 1. Copier le dossier "templates"
    src_templates = joinpath(src_cluster_path, "templates")
    dest_templates = joinpath(dest_cluster_path, "templates")
    mkpath(dest_templates)
    for pdb_file in readdir(src_templates)
        cp(joinpath(src_templates, pdb_file), joinpath(dest_templates, pdb_file); force=true)
    end
    # 2. Appliquer MAFFT sur sequences.a3m
    src_seq_file = joinpath(src_cluster_path, "sequences.a3m")
    dest_seq_file = joinpath(dest_cluster_path, "sequences.a3m")

    if isfile(src_seq_file)
        run(pipeline(`$(MAFFT_jll.mafft()) --quiet $src_seq_file`, stdout=dest_seq_file))
    else
        @warn "Fichier sequences.a3m manquant dans $src_cluster_path"
    end
    msa = read(dest_seq_file, A3M)
    output_file=adjustreference!(msa)
    write(dest_seq_file,output_file,A3M)
end

AlphaConformers.run_alphafold(dest_root, colabfold_path=COLABFOLD_PATH)
