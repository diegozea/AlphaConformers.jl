#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

#=
#PBS -l host=node48
#PBS -l walltime=900:00:00
#PBS -l mem=100gb
#PBS -l ncpus=40
#PBS -j oe
=#

#=
using Revise
using CSV, DataFrames
using AlphaConformers

const PATH = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/mohit"
cd(PATH)

const PDB_FOLDER = nothing

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"

output_dir = joinpath(PATH, "first_run")
if isdir(output_dir)
    rm(output_dir; recursive=true, force=true)
end
mkdir(output_dir)

const REF_PDB = joinpath(PATH, "P51993_35_359.pdb")

AlphaConformers.alphaconformers(REF_PDB, FOLDSEEK_DB, ALPHAFOLD_DB, PDB_FOLDER, output_dir)
AlphaConformers.run_alphafold(output_dir, colabfold_path=COLABFOLD_PATH)
=#

using Revise
using CSV, DataFrames
using AlphaConformers

# Define constants and paths
const PATH = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/mohit"
cd(PATH)

const REF_PDB = joinpath(PATH, "P51993_35_359.pdb")
const PDB_FOLDER = nothing  # Adjust if you have a specific PDB folder
const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"
const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"
const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"
const EVALUE_CUTOFF = 1e-5  # Adjust the e-value cutoff if needed

# Set up the output directory
output_dir = joinpath(PATH, "first_run")
if isdir(output_dir)
    rm(output_dir; recursive=true, force=true)
end
mkdir(output_dir)

# Step 1: Run Foldseek
@info "Running Foldseek"
foldseek_output = AlphaConformers.run_foldseek(
    REF_PDB,
    join([FOLDSEEK_DB, ALPHAFOLD_DB], ","),
    out_folder=output_dir
)
println("Foldseek output:", foldseek_output)

# Step 2: Merge Foldseek tables
@info "Merging Foldseek tables"
merged_table = AlphaConformers.merge_tables([
    foldseek_output[1].table_file,
    foldseek_output[2].table_file
])
println("Merged table before filtering:", merged_table)

# Step 3: Filter the merged table by e-value cutoff
@info "Filtering merged table by e-value"
filter!(row -> row.evalue < EVALUE_CUTOFF, merged_table)
println("Merged table after filtering:", merged_table)

# Step 4: Merge MSAs
@info "Merging MSAs"
merged_msa = AlphaConformers.merge_msas(merged_table)
println("Merged MSA:", merged_msa)

# Step 5: Get the aligned structures
@info "Getting the aligned structures"
structures = AlphaConformers._get_aligned_structures(merged_table)

# Step 6: Add known conformations
@info "Adding known conformations"
sifts_uniprot_mapping = AlphaConformers.get_uniprot_mapping()
target2uniprot, expanded_table = AlphaConformers.add_known_conformations!(
    deepcopy(merged_table),
    sifts_uniprot_mapping
)

# Step 7: Get uniprot to targets mapping
uniprot2targets = AlphaConformers.get_uniprot2targets(target2uniprot, expanded_table)

# Step 8: Measure RMSDs
@info "Measuring RMSDs"
rmsds = AlphaConformers.process_known_conformations!(
    structures,
    expanded_table,
    target2uniprot,
    uniprot2targets,
    pdb_folder=PDB_FOLDER
)

# Fill RMSDs in the expanded table
AlphaConformers.fill_rmsds!(rmsds, expanded_table, structures)

# Step 9: Cluster structures
@info "Clustering structures"
large_small_pairs, large_cl2msa, small_cl2pdb = AlphaConformers.create_template_clusters(
    rmsds,
    expanded_table,
    merged_msa,
    structures
)

# Step 10: Create folder structure
@info "Creating folder structure"
AlphaConformers.create_folder_structure(
    large_small_pairs,
    large_cl2msa,
    small_cl2pdb,
    out_folder=output_dir
)
