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

# Load necessary packages on all workers
using DataFrames, CSV
import MIToS
using MIToS.MSA
using MIToS.PDB
using MMseqs2_jll


############################################################ Functions ##############################################################


function write_a3m(filename, fasta)
    open(filename, "w") do io
        for (header, seq) in fasta
            println(io, ">$header")
            println(io, seq)
        end
    end
end
####################################################### MAIN ##########################################################

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

# 1. Charger le CSV
@info "Start"
model = "1"
file_path_df_final="pdb_information_details_final_mutation_cluster_reformatted_filter_foldseek2_final2_1.csv"
df=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND

# 2. Filtrer les lignes où LIGANDS est manquant
df_holo = filter(row -> ismissing(row.LIGANDS), df)

fasta = []
for row in eachrow(df_holo)
    holo_id = String(row.UNIPROT)
    path_file=joinpath(PDB_FOLDER,uppercase(String(row.PDB))*".pdb")
    if !isfile(path_file)
        MIToS.PDB.downloadpdb(String(row.PDB); format=MIToS.PDB.PDBFile, filename=path_file)
    end
    try 
        pdb_file = MIToS.PDB.read_file(path_file, MIToS.PDB.PDBFile)
        if pdb_file === nothing
            @warn "PDB file is empty or could not be read: $path_file"
            continue
        end

        chain = String(row.CHAIN)
        sequences=MIToS.PDB.modelled_sequences(pdb_file;chain=chain)
        @show sequences
        seq = sequences[(model=model, chain=chain)]        
        @show seq
        push!(fasta, (holo_id*"_"*chain, seq))
    
    catch e
        @warn "Failed to read PDB file: $path_file, error: $(e)"
        continue
    end
    
end

# 3. Écrire le fichier FASTA

write_a3m("dataset_holo_MMSEQ.fasta", fasta)

run(`$(MMseqs2_jll.mmseqs()) easy-cluster dataset_holo_MMSEQ.fasta homologs tmp --min-seq-id 0.3 -c 0.8 --alignment-mode 3 --cluster-mode 1`)

clusters = CSV.read("homologs_cluster.tsv", DataFrame, header=false, delim='\t')
@show first(clusters,10)
@info "End !"
########################################################## END ##########################################################