# Pipeline
# ========
# 1. Take a protein sequence and one of its conformations
#    In the future, I can think of doing a first AlphaFold model to use as the starting point.
# 2. Look for proteins sharing that conformation using the FoldSeek API.
# 3. For each of those proteins, look alternative conformations or structural models.
#    For that, I can use all the structures from protein with a very high 
#    percentage of sequence identity (that should be a parameter set a 99% by default).
# 4. Cluster all the structures, and define set of possible conformations for the protein.
# 5. Use the AlphaFold model to predict the structure of the protein using the set of 
#    conformations as templates.
# 6. Compare the predicted structure with the other structures to assess wheter the desired
#    conformation was modeled correctly.

# qsub -I -l host=node48 -l walltime=300:00:00 -l ncpus=20 -l mem=64g

using Pkg
# Pkg.activate(abspath(@__DIR__, ".."))
Pkg.activate("/home/diego.zea/.julia/dev/AlphaConformers/")

using MIToS.PDB
using CSV
using DataFrames

const WORKING_DIR = "/alpha/runs/diego.zea/first_test"

cd(WORKING_DIR)

pdb_code, chain_code = "1EX6", "B"

# pdb_file = downloadpdb(pdb_code, format=PDBFile)
chain = read("$(pdb_code).pdb.gz", PDBFile, chain=chain_code)
pdb_file = abspath("$(pdb_code)_$(chain_code).pdb")
write("$(pdb_code)_$(chain_code).pdb", chain, PDBFile)

# run foldseek locally

function run_foldseek_search(pdb_file)
    foldseek = "/store/EQUIPES/AMIG/PROGRAMMES/foldseek/bin/foldseek"
    # db = "/store/EQUIPES/AMIG/MEMBERS/carla.martins/foldseek/db_pdb/pdb"
    db = "/alpha/database/pdb"
    tmp_folder = mktempdir()
    out_file = "$(pdb_file)_results.m8"
    run(`$foldseek easy-search $pdb_file $db $out_file $tmp_folder`)
    out_file
end

search_results_file = run_foldseek_search(pdb_file)

# parse foldseek results

column_names = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits"]

search_results = DataFrame(CSV.File(search_results_file, delim='\t', header=column_names))

# ---

# do the evalue allow to cluster conformations?
import UnicodePlots

UnicodePlots.scatterplot(search_results.fident, log10.(search_results.evalue), canvas=UnicodePlots.AsciiCanvas, border=:ascii)

# ---

# cluster structures using foldseek and see if the conformations are located in different clusters

function run_foldseek_clustering(search_results::DataFrame)
    # PATH to the PDB files
    pdb_folder = "/alpha/database/pdb/pdb_files"

    # Create a local folder for PDB files
    if isdir("first_targets")
        rm("first_targets"; recursive=true)
    end
    mkdir("first_targets")
    local_pdb_folder = abspath("first_targets")
    
    # Go through the target column of search results
    for row in eachrow(search_results)
        target = row[:target]
        # Split target into pdb_code and chain_code
        fields = split(target, '_')
        pdb_code = fields[1]
        if length(fields) == 2
            chain_code = fields[2]
        else
            chain_code = nothing
        end
        # Define original and local file paths
        original_pdb_path = joinpath(pdb_folder, pdb_code)
        local_pdb_path = joinpath(local_pdb_folder, pdb_code)
        if chain_code !== nothing
            # Read only specific chain and write to local file
            chain = read(original_pdb_path, PDBFile, chain=String(chain_code))
            write(local_pdb_path, chain, PDBFile)
        else
            # Create a symlink to the original file
            symlink(original_pdb_path, local_pdb_path)
        end
    end
    
    # Run foldseek easy-cluster
    tmp_folder = mktempdir()
    out_file = "first_targets"
    foldseek = "/store/EQUIPES/AMIG/PROGRAMMES/foldseek/bin/foldseek"
    run(`$foldseek easy-cluster $local_pdb_folder $out_file $tmp_folder --alignment-mode 3 --tmalign-fast 1 --rescore-mode 4 --cov-mode 1 --cluster-mode 1 --similarity-type 1 --min-seq-id 0.0 -e 0.0000001`)
    return out_file
end

run_foldseek_clustering(search_results)

# I made just a few attempts with one example, but it seems that foldseek cannot 
# effectively distinguish between conformations belonging to different clusters.

#---

# delete the query if it is in the target column

function delete_query_from_target!(search_results)
    pdb = "$(pdb_code).pdb"
    pdb_chain = "$(pdb_code).pdb_$(chain_code)"
    filter!(row -> row.target != pdb && row.target != pdb_chain, search_results)
end

delete_query_from_target!(search_results)

# delete all the rows with a sequence identity higher than 98%
filter!(row -> row.fident < 0.98, search_results)

# TODO: Decide whether to use UniProt accessions to remove known conformations in the 
# final pipeline. Consider making it optional.

