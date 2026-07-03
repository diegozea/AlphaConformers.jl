module AlphaConformers

import BioAlignments
import BioStructures
import MIToS
import Foldseek_jll
import MAFFT_jll
import USalign_jll
import CSV
import DataFrames
import Clustering
import Downloads
import OrderedCollections
import PairwiseListMatrices
import Distributed
import ProgressMeter
import Combinatorics
import Random
import Statistics
import MultivariateStats
import Distances
import Plots as plt
import StatsPlots
import HTTP
import JSON3
import YAML
using Glob: glob
using DataFrames: DataFrame, nrow



export foldseek_search, # foldseek.jl
    read_foldseek_search_results,
    run_foldseek,
    merge_tables,
    merge_msas,
    get_aligned_positions,
    create_pdb_folder, # pdb_folders.jl
    usalign, # usalign.jl
    get_uniprot_mapping, # sifts.jl
    get_pfam_mapping,
    get_uniprot_acc,
    get_pdb_codes,
    delete_query_from_target,
    get_residues_and_sequence, # af_input.jl
    save_sequences,
    align_sequences,
    clean_msa,
    create_pdb_lists,
    create_msa_and_templates,
    organize_files,
    run_alphafold3,
    organize_files_af3,
    run_boltz2,
    organize_files_boltz,
    read_a3m,
    clean_msa_template_names,
    structural_alignment,
    _read_pdb_chain,
    get_msa_sequence_afdb,
    run_alphafold_input_structure,
    found_best_prediction,
    triage_outputs,
    alphaconformers,
    PreparedInputs,
    prepare_inputs,
    add_known_conformations!,
    create_template_clusters_hobohm,
    create_folder_structure_hobohm,
    run_cmd,
    run_alphafold,
    run_deepaccnet,
    cluster_conformers

include("utils.jl")
include("foldseek.jl")
include("clustering.jl")
include("sifts.jl")
include("pdb_folders.jl")
include("usalign.jl")
include("af_input.jl")
include("run_alphafold.jl")
include("organize_files.jl")
include("run_alphaconformers.jl")
include("run_deepaccnet.jl")
include("conformer_clustering.jl")
include("conformer_clustering_plots.jl")

end
