module AlphaConformers

import MIToS
import Foldseek_jll
import MAFFT_jll
import USalign_jll
import CSV
import DataFrames
import Clustering
import Downloads

using TestItems

export  foldseek_search, # foldseek.jl
        read_foldseek_search_results,
        create_pdb_folder, # pdb_folders.jl
        usalign, # usalign.jl
        structural_clustering, # clustering.jl
        get_uniprot_mapping, # sifts.jl
        get_uniprot_acc,
        get_pdb_codes,
        delete_query_from_target!,
        list_known_conformations,
        get_residues_and_sequence, # af_input.jl
        save_sequences,
        align_sequences,
        clean_msa,
        create_pdb_lists,
        create_alpha_fold_inputs


include("foldseek.jl")
include("pdb_folders.jl")
include("usalign.jl")
include("clustering.jl")
include("sifts.jl")
include("af_input.jl")

end
