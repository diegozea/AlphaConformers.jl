#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10
#SBATCH --output=train_valid_test-%j.out
=#
import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")

using DataFrames, CSV
import MIToS
using MIToS.MSA
using MIToS.PDB
using MMseqs2_jll
const PDB_FOLDER = "/alpha/database/pdb/pdb_files"
# 1. Charger le CSV
file_path_df_final="pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek2_final2_1.csv"
df=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND

# 2. Filtrer les lignes où LIGANDS est manquant
df_holo = filter(row -> ismissing(row.LIGANDS), df)

fasta = String[]
for row in eachrow(df_holo)
    holo_id = row.UNIPROT
    pdb_file = joinpath(PDB_FOLDER, String(row.PDB)*".pdb.gz")
    if !isfile(pdb_file)
        MIToS.PDB.downloadpdb(String(row.PDB); format=MIToS.PDB.PDBFile, filename=pdb_file)
    end
    if pdb_file === nothing
        continue
    end
    structure = read(pdb_file, PDBFile)
    chain = structure[String(row.CHAIN)]
    sequence = MIToS.PDB.sequence(chain)
    push!(fasta, ">$(holo_id)\n$(sequence)")
end
# 3. Écrire le fichier FASTA
msa_file = "uniprot_sequences.fasta"
MIToS.MSA.write_file(msa_file, fasta, MIToS.MSA.FASTA)

run(`$(MMseqs2_jll.mmseqs()) easy-cluster uniprot_sequences.fasta homologs tmp --min-seq-id 0.3 -c 0.8 --alignment-mode 3 --cluster-mode 1`)

clusters = CSV.read("clusters.tsv", DataFrame, header=false, delim='\t')
@show first(clusters,10)