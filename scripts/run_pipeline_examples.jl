#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia110

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

const PATH = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/talk_subset"
cd(PATH)

const DATA = DataFrame(CSV.File(joinpath(PATH, "selected_examples.csv")))

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const ALPHAFOLD_DB = "/alpha/database/afdb/afdb_up"

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"


# ┌ Error: Error processing 1AKZ_A: DimensionMismatch("arrays could not be broadcast to a common size; got a dimension with lengths 122 and 117")
# ┌ Error: Error processing 1HW1_B: ErrorException("Inconsistent dictionary sizes")
# ┌ Error: Error processing 1GQN_A: HTTP.Exceptions.StatusError(404, "GET", "/download/7QXX.pdb.gz", HTTP.Messages.Response:
# ┌ Error: Error processing 2LHS-6_A: HTTP.Exceptions.StatusError(404, "GET", "/download/7SOE.pdb.gz", HTTP.Messages.Response:
# ┌ Error: Error processing 1FMF-4_A: DimensionMismatch("arrays could not be broadcast to a common size; got a dimension with lengths 151 and 150")
# ┌ Error: Error processing 2F63-4_A: ErrorException("Inconsistent dictionary sizes")
# ┌ Error: Error processing 2UZ5-9_A: ErrorException("Inconsistent dictionary sizes")
# ┌ Error: Error processing 4AKE_B: DimensionMismatch("arrays could not be broadcast to a common size; got a dimension with lengths 42 and 41")
# ┌ Error: Error processing 1O1U-9_A: ErrorException("Inconsistent dictionary sizes")
# ┌ Error: Error processing 6OY9_B: HTTP.Exceptions.StatusError(404, "GET", "/download/4U3N.pdb.gz", HTTP.Messages.Response:
# ┌ Error: Error processing 1MUT-11_A: HTTP.Exceptions.StatusError(404, "GET", "/download/7Q93.pdb.gz", HTTP.Messages.Response:
# ┌ Error: Error processing 4LP5_A: ErrorException("Inconsistent dictionary sizes")

#=

row = DATA[1, :]

@info "Processing $(row.apo_id)"   
apo_id = row.apo_id
apo_dir = joinpath(PATH, apo_id)
if isdir(apo_dir)
    rm(apo_dir; recursive=true, force=true)
end
mkdir(apo_dir)

m = match(r"([A-Za-z0-9]{4})-?([0-9]*)_([A-Z0-9a-z]+)", apo_id)
ref_pdb_code = String(m.captures[1])
ref_chain = String(m.captures[3])
ref_model = isempty(m.captures[2]) ? "1" : String(m.captures[2])
ref_pdb = joinpath(PDB_FOLDER, "$ref_pdb_code.pdb")

AlphaConformers.alphaconformers(ref_pdb, FOLDSEEK_DB, ALPHAFOLD_DB, PDB_FOLDER, apo_dir)
AlphaConformers.run_alphafold(apo_dir, colabfold_path=COLABFOLD_PATH)


=#

for row in eachrow(DATA)
    try
        @info "Processing $(row.apo_id)"   
        apo_id = row.apo_id
        apo_dir = joinpath(PATH, apo_id)
        if isdir(apo_dir)
            rm(apo_dir; recursive=true, force=true)
        end
        mkdir(apo_dir)

        m = match(r"([A-Za-z0-9]{4})-?([0-9]*)_([A-Z0-9a-z]+)", apo_id)
        ref_pdb_code = String(m.captures[1])
        ref_chain = String(m.captures[3])
        ref_model = isempty(m.captures[2]) ? "1" : String(m.captures[2])
        ref_pdb = joinpath(PDB_FOLDER, "$ref_pdb_code.pdb")
        
        AlphaConformers.alphaconformers(ref_pdb, FOLDSEEK_DB, ALPHAFOLD_DB, PDB_FOLDER, apo_dir)
        AlphaConformers.run_alphafold(apo_dir, colabfold_path=COLABFOLD_PATH)

        # paths = create_alpha_fold_inputs(apo_dir, ref_pdb, ref_chain, ref_model,
        #     foldseek_db=FOLDSEEK_DB, pdb_db=PDB_FOLDER, testing=true)
        # run_alphafold(paths.clusters, colabfold_path=COLABFOLD_PATH)
    catch err
        @error "Error processing $(row.apo_id): $(err)"
    end
end


#= 
ERROR: InexactError: Int64(10.63014581273465)
Stacktrace:
  [1] Int64
    @ ./float.jl:909 [inlined]
  [2] _nelements
    @ ~/.julia/packages/PairwiseListMatrices/OFqua/src/pairwiselistmatrix.jl:278 [inlined]
  [3] PairwiseListMatrices.PairwiseListMatrix(list::Vector{Float64}, diagonal::Bool, diagonalvalue::Float64)
    @ PairwiseListMatrices ~/.julia/packages/PairwiseListMatrices/OFqua/src/pairwiselistmatrix.jl:296
  [4] from_table(table::DataFrame, diagonal::Bool; labelcols::Vector{Int64}, valuecol::Int64, diagonalvalue::Float64)
    @ PairwiseListMatrices ~/.julia/packages/PairwiseListMatrices/OFqua/src/pairwiselistmatrix.jl:1179
  [5] from_table
    @ ~/.julia/packages/PairwiseListMatrices/OFqua/src/pairwiselistmatrix.jl:1169 [inlined]
  [6] get_rmsd_matrix(rmsds::Dict{Tuple{String, String}, Float64}, targets::Set{String})
    @ AlphaConformers ~/.julia/dev/AlphaConformers/src/foldseek.jl:507
  [7] cluster_structures(rmsds::Dict{Tuple{String, String}, Float64}, targets::Set{String})
    @ AlphaConformers ~/.julia/dev/AlphaConformers/src/foldseek.jl:515
  [8] create_template_clusters(rmsds::Dict{…}, expanded_table::DataFrame, msa::MIToS.MSA.AnnotatedMultipleSequenceAlignment, structures::OrderedCollections.OrderedDict{…})
    @ AlphaConformers ~/.julia/dev/AlphaConformers/src/foldseek.jl:579
  [9] alphaconformers(input_pdb::String, pdb_db::String, alphafold_db::String, pdb_folder::String, out_folder::String; evalue_cutoff::Float64)
    @ AlphaConformers ~/.julia/dev/AlphaConformers/src/foldseek.jl:634
 [10] alphaconformers(input_pdb::String, pdb_db::String, alphafold_db::String, pdb_folder::String, out_folder::String)
    @ AlphaConformers ~/.julia/dev/AlphaConformers/src/foldseek.jl:612
 [11] top-level scope
    @ REPL[882]:1
Some type information was truncated. Use `show(err)` to see complete types.
=#