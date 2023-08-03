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

ENV["PATH"] *= ":/usr/local/bin"
ENV["SIFTS_DB"] = "/home/diego.zea"

const _PATH = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset"
cd(_PATH)

if !isdir("tmp")
    mkdir("tmp")
end

cd("tmp")

const PATH = pwd()

const DATA = DataFrame(CSV.File(joinpath(_PATH, "selected_examples.csv")))

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"

# apo_id = "6OY9_B"
# apo_id = "2LHS-6_A"

function run_alphaconformers(apo_id::String)
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
end

# const RERUN = ["1FMF-4_A", "2F63-4_A", "1O1U-9_A", "4LP5_A"] 

# ERROR: AssertionError: File /store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/tmp/1O1U-9_A/pdb/7FY8.pdb_F does not exist
# Stacktrace:
#   [1] usalign(pdb_file_a::String, pdb_folder::String, pdb_list::Vector{String}; kargs::Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}})
#     @ AlphaConformers ~/.julia/dev/AlphaConformers/src/usalign.jl:104
#   [2] usalign
#     @ ~/.julia/dev/AlphaConformers/src/usalign.jl:103 [inlined]
#   [3] structural_clustering(query_pdb::String, pdb_folder::String, targets::Vector{String}; rmsd_cutoff::Float64)
#     @ AlphaConformers ~/.julia/dev/AlphaConformers/src/clustering.jl:66

# map(run_alphaconformers, RERUN)
