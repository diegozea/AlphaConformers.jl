#!/home/julie.daniel/.julia/juliaup/julia-1.11.7+0.x64.linux.gnu/bin/julia

#=
#SBATCH --nodelist=node48
#SBATCH --time=900:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --output=run_AlphaConformers.jl.o%j.out
=#

import Pkg
Pkg.activate(joinpath(@__DIR__, "update_MIToS_321"))

using AlphaConformers
using MIToS
using MIToS.PDB
using DataFrames
import CSV
using Revise
using Glob
using Statistics
using ArgParse

################################################## FUNCTIONS ######################################################

function parse_commandline()
    s = ArgParseSettings(
        description = "Run AlphaConformers pipeline on a PDB file: structural search, clustering, and AlphaFold prediction.",
    )

    @add_arg_table! s begin
        "--apo_pdb", "-q"
        help = "PDB ID of the apo form to use as query (e.g. 6r17)"
        arg_type = String
        required = true
        "--apo_chain", "-c"
        help = "Chain ID of the apo form (e.g. C)"
        arg_type = String
        required = true
        "--path", "-p"
        help = "Path to the main directory containing the apo/holo files and where results will be saved"
        arg_type = String
        required = true
        "--n_threads", "-t"
        help = "Number of threads to use for AlphaConformers (default: all available)"
        arg_type = Int
        default = Threads.nthreads()
        "--pdb_folder"
        help = "Path to the PDB mmCIF files"
        arg_type = String
        default = "/alpha/database/pdb/mmcif_files"
        "--foldseek_db"
        help = "Path to the Foldseek database"
        arg_type = String
        default = "/alpha/database/pdb/fullpdb_mmcif_files"
        "--alphafold_db"
        help = "Path to the AlphaFold database"
        arg_type = String
        default = "/alpha/database/afdb_v6/fullafdb_v6"
        "--sif_path"
        help = "Path to the ColabFold Singularity image"
        arg_type = String
        default = expanduser(
            "/store/EQUIPES/AMIG/SCRIPTS/sif_images/ColabFold_AF2_1-5-5/colabfold-1.5.5-cuda12.2.2.sif",
        )
        "--cache_dir"
        help = "Path to the ColabFold cache directory"
        arg_type = String
        default = "/store/EQUIPES/AMIG/SCRIPTS/sif_images/ColabFold_AF2_1-5-5/cache"
        "--evalue_cutoff"
        help = "E-value cutoff for Foldseek (default: NaN = no cutoff)"
        arg_type = Float64
        default = NaN
        "--cutoff"
        help = "Clustering cutoff (default: 1.0)"
        arg_type = Float64
        default = 1.0
        "--overwrite"
        help = "Overwrite the output directory if it already exists. Without this flag, the script will error if the folder exists."
        action = :store_true
        "--test_analyse"
        help = "Delete structure found by Foldseek that have the same UNIPROT ID as the query, to test the analyse without structure from the same proteins"
        action = :store_false
        "--keep_query"
        help = "Keep the query structure in the Foldseek results (default: true). Set to false to remove it and test the analysis without the query structure. Available only if --test_analyse is set to true."
        action = :store_true
        "--full_seq"
        help = "Specify a full sequence to use instead of the query .pdb sequence. The sequence should be provided as a string of amino acid one-letter codes (e.g. 'MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFP'). If missing, the script will use the sequence from the query .pdb file. This option is useful if the query .pdb file is incomplete or if you want to test the analysis with a different sequence."
        arg_type = String
        default = ""

    end

    return parse_args(s)
end

"""
    run_alphaconformers(apo_pdb, apo_chain, path, n_threads, pdb_folder, foldseek_db,
                        alphafold_db, sif_path, cache_dir; evalue_cutoff, cutoff)

Execute the AlphaConformers pipeline on a given apo PDB structure.

Runs structural search with Foldseek, clusters the results, then launches
AlphaFold predictions for each cluster via ColabFold.

# Arguments
- `apo_pdb::String`       : PDB ID of the apo form (e.g. `"6r17"`).
- `apo_chain::String`     : Chain ID of the apo form (e.g. `"C"`).
- `path::String`          : Main directory for input files and results output.
- `n_threads::Int`        : Number of threads to use.
- `pdb_folder::String`    : Path to the PDB mmCIF files.
- `foldseek_db::String`   : Path to the Foldseek database.
- `alphafold_db::String`  : Path to the AlphaFold database.
- `sif_path::String`      : Path to the ColabFold Singularity image.
- `cache_dir::String`     : Path to the ColabFold cache directory.
- `evalue_cutoff::Float64`: E-value cutoff for Foldseek (`NaN` = no cutoff).
- `cutoff::Float64`       : Clustering cutoff value (default `1.0`).

# Output
- One output folder per PDB file with AlphaConformers results.
- Each cluster subfolder contains AlphaFold prediction results.

# Notes
- AlphaConformers takes ~40 min per PDB, then waits for AlphaFold per cluster.
- Total runtime depends on the number of clusters found.
- To reduce cluster count, tune `evalue_cutoff` and `cutoff`.
"""
function run_alphaconformers(
    apo_pdb::String,
    apo_chain::String,
    path::String,
    n_threads::Int,
    pdb_folder::String,
    foldseek_db::String,
    alphafold_db::String,
    sif_path::String,
    cache_dir::String;
    evalue_cutoff::Float64 = NaN,
    cutoff::Float64 = 1.0,
    overwrite::Bool = false,
    test_analyse::Bool = false,
    keep_query::Bool = true,
    full_seq::Union{String,Missing} = missing,
)
    cd(path)

    db = [foldseek_db, alphafold_db]

    filename = string(apo_pdb, "_", apo_chain, ".pdb")
    ref_pdb = joinpath(path, filename)
    @show ref_pdb

    output_dir = joinpath(path, apo_pdb * "_AlphaConformer")

    if isdir(output_dir)
        if !overwrite
            error("""
            Output directory already exists: $output_dir
            Re-run with --overwrite to delete it and start fresh.
            Warning: all existing data in this folder will be permanently deleted.
            """)
        end
        @warn "Overwriting existing directory: $output_dir"
        rm(output_dir; recursive = true, force = true)
    end
    mkdir(output_dir)
    println(output_dir)

    AlphaConformers.alphaconformers(
        ref_pdb,
        pdb_folder,
        output_dir;
        n_threads = n_threads,
        db = db,
        evalue_cutoff = evalue_cutoff,
        cutoff = cutoff,
        test_analyse = false,
        keep_query = true,
        full_seq = missing,
    )
    AlphaConformers.run_alphafold_one_run(output_dir, sif_path, cache_dir)

    @info "Done"
end

################################################## MAIN ######################################################
#Run the pipeline with the command line arguments
#julia run_AlphaConformers.jl \
#    --apo_pdb 6r17 \
#    --apo_chain C \
#    --path /store/.../data/ \
#   --foldseek_db /autre/chemin/foldseek \
#    --cutoff 0.8 \
#    --evalue_cutoff 0.001
#    --overwrite

args = parse_commandline()

run_alphaconformers(
    args["apo_pdb"],
    args["apo_chain"],
    args["path"],
    args["n_threads"],
    args["pdb_folder"],
    args["foldseek_db"],
    args["alphafold_db"],
    args["sif_path"],
    args["cache_dir"];
    evalue_cutoff = args["evalue_cutoff"],
    cutoff = args["cutoff"],
    overwrite = args["overwrite"],
    test_analyse = args["test_analyse"],
    keep_query = args["keep_query"],
    full_seq = args["full_seq"],
)
