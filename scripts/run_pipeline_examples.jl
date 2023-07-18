#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l host=node48
#PBS -l walltime=900:00:00
#PBS -l mem=100gb
#PBS -l ncpus=40
#PBS -j oe
=#

using Revise
using CSV, DataFrames
using AlphaConformers

const PATH = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset"
cd(PATH)

const DATA = DataFrame(CSV.File(joinpath(PATH, "selected_examples.csv")))

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"

for row in eachrow(DATA)
    try
        @info "Processing $(row.apo_id)"   
        apo_id = row.apo_id
        apo_dir = joinpath(PATH, apo_id)
        isdir(apo_dir) || mkdir(apo_dir)

        m = match(r"([A-Za-z0-9]{4})-?([0-9]*)_([A-Z0-9a-z]+)", apo_id)
        ref_pdb_code = String(m.captures[1])
        ref_chain = String(m.captures[3])
        ref_model = isempty(m.captures[2]) ? "1" : String(m.captures[2])
        ref_pdb = joinpath(PDB_FOLDER, "$ref_pdb_code.pdb")

        paths = create_alpha_fold_inputs(apo_dir, ref_pdb, ref_chain, ref_model,
            foldseek_db=FOLDSEEK_DB, pdb_db=PDB_FOLDER, testing=true)

        run_alphafold(paths.clusters, colabfold_path=COLABFOLD_PATH)
    catch err
        @error "Error processing $(row.apo_id): $(err)"
    end
end
