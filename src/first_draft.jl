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

using MIToS.PDB
using CSV
using DataFrames


pdb_code, chain_code = "1EX6", "B"

# pdb_file = downloadpdb(pdb_code, format=PDBFile)
chain = read("test/data/$(pdb_code).pdb.gz", PDBFile, chain=chain_code)
pdb_file = abspath("test/data/$(pdb_code)_$(chain_code).pdb")
write("test/data/$(pdb_code)_$(chain_code).pdb", chain, PDBFile)

# run foldseek locally

function run_foldseek_search(pdb_file)
    foldseek = "/store/EQUIPES/AMIG/PROGRAMMES/foldseek/bin/foldseek"
    db = "/store/EQUIPES/AMIG/MEMBERS/carla.martins/foldseek/db_pdb/pdb"
    tmp_folder = mktempdir()
    out_file = "$(pdb_file)_results.m8"
    run(`$foldseek easy-search $pdb_file $db $out_file $tmp_folder`)
    out_file
end

search_results_file = run_foldseek_search(pdb_file)

# parse foldseek results

column_names = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits"]

search_results = DataFrame(CSV.File(search_results_file, delim='\t', header=column_names))

# delete the query if it is in the target column
query = lowercase("$(pdb_code)_$(chain_code)")
filter!(row -> row.target != query, search_results)

