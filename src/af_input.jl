# These functions create the folder structure defining AlphaFold 2 inputs, including the 
# the MSA, and the templates.

"""
    get_residues_and_sequence(pdb_file; chain=MIToS.PDB.All, model="1")

This function reads a pdb file and returns the residues and the sequence of the chain
specified by the `chain` argument. The `chain` keyword argument can be a string with the 
chain name or the `MIToS.PDB.All` type. However, this function will throw an error if more 
than one chain is selected; so `All` is used only to avoid specifiying a chain when there is 
only one in the file. The `model` argument is used to select the model to read (by default
the first model is selected).
"""
function get_residues_and_sequence(pdb_file; 
        chain::Union{Type{MIToS.PDB.All}, String}=MIToS.PDB.All,
        model::String="1")
    query_res = read(pdb_file, MIToS.PDB.PDBFile, 
        model=model, chain=chain, group="ATOM",
        onlyheavy=true, occupancyfilter=true)
    sequences = collect(values(MIToS.PDB.modelled_sequences(query_res)))
    @assert length(sequences) == 1 "The query must be a single chain"
    (residues=query_res, sequence=sequences[1])
end


"""
    save_sequences(cluster_folder::String, pdb_files::Vector{String}; chains, models)

This function uses the `get_residues_and_sequence` function to create a fasta file with 
the sequences of the structures in the cluster. The second arguments should be a vector
with the paths to the pdb files of the structures in the cluster. The `chains` and `models`
arguments are vectors with the chain and model to read for each pdb file. By default, the
first model and all the chains are selected.
"""
function save_sequences(cluster_folder::String, pdb_files::Vector{String}; 
        chains=fill(MIToS.PDB.All, length(pdb_files)),
        models::Vector{String}=fill("1", length(pdb_files)))
    open(joinpath(cluster_folder, "sequences.fasta"), "w") do file 
        for (pdb_file, chain, model) in zip(pdb_files, chains, models)
            pdb_name = basename(pdb_file)
            pdb_data = get_residues_and_sequence(pdb_file; chain=chain, model=model)
            println(file, ">$pdb_name")
            println(file, pdb_data.sequence)
        end
    end
end

"""
    align_sequences(sequences::String)

This function uses `MAFFT` (thanks to the `MAFFT_jll` package) to align the sequences in 
the given fasta file. It returns a MSA object from `MIToS.MSA` that keeps the original 
order of the sequences.
"""
function align_sequences(sequences::String)
    @assert isfile(sequences) "$(sequences) does not exist."
    mktemp() do aln_path, _
        run(pipeline(`$(MAFFT_jll.mafft()) --quiet $sequences`, stdout=aln_path))
        @assert isfile(aln_path) && filesize(aln_path) > 0 "The MSA was not properly created."
        read(aln_path, MIToS.MSA.FASTA)
    end
end


"""
    clean_msa(msa::MIToS.MSA.AnnotatedMultipleSequenceAlignment)

This function removes the sequences with more than 50% of gaps and the columns with gaps 
in the reference sequence (the first one). The function returns a new MSA object.
"""
function clean_msa(msa::MIToS.MSA.AnnotatedMultipleSequenceAlignment)
    msa_ref = MIToS.MSA.adjustreference(msa)
    msa_ref[vec(MIToS.MSA.coverage(msa_ref) .≥ 0.5), :]
end


# NOTE: MODELS: We want to select a specific model for the query, as we usually consider
# all the models when looking for the pair of conformations with the highest RMSD 
# (maximum conformational diversity as in CoDNaS). However, we will use only the first 
# model for the templates, to reduce the number of structural alignments. This could be
# changed in the future, but it is not a priority now.


"""
    create_pdb_lists(ref_pdb, ref_chain, ref_model, pdb_files, chains, models)

This function creates the lists of pdb files, chains and models to use as input for
AlphaFold 2; defining the sequence and templates to consider.
"""
function create_pdb_lists(ref_pdb, ref_chain, ref_model, pdb_files, chains, models)
    ref_abspath = abspath(ref_pdb)
    pdb_abspaths = abspath.(pdb_files)
    chains = deepcopy(chains)
    models = deepcopy(models)
    # Check that the query is the first pdb file in the list
    if abspath(first(pdb_abspaths)) != ref_abspath
        # if not, check that it is not in another position
        ref_pos = findfirst(==(ref_abspath), pdb_abspaths)
        if isnothing(ref_pos)
            # if missing, add it as the first element of the list
            pushfirst!(pdb_abspaths, ref_abspath)
            pushfirst!(chains, ref_chain)
            pushfirst!(models, ref_model)
        else
            # otherwise, move it to the first position and inform the user with a warning
            pdb_abspaths[1], pdb_abspaths[ref_pos] = pdb_abspaths[ref_pos], pdb_abspaths[1]
            chains[1], chains[ref_pos] = chains[ref_pos], chains[1]
            models[1], models[ref_pos] = models[ref_pos], models[1]
        end
    end
    # Check that the number of pdb files, chains and models is the same after the changes
    @assert length(pdb_abspaths) == length(chains) == length(models) "The number of pdb files, chains and models must be the same"

    (;pdb_files=pdb_abspaths, chains, models)
end

"""
    create_alpha_fold_inputs(cluster_folder, ref_pdb, ref_chain, ref_model, 
        pdb_files, chains, models)

This function creates the folder structure defining AlphaFold 2 inputs, including the
MSA and the templates. The `cluster_folder` argument is the path to the folder where the
inputs will be created. The `ref_pdb`, `ref_chain` and `ref_model` arguments define the
query structure. The `pdb_files`, `chains` and `models` arguments define the structures
to use as templates. Note that the query structure or reference can also be included in 
`pdb_files`, `chains` and `models`; in that case, the `ref_chain` and `ref_model` arguments
are ignored. The function returns the path to the pdb files, and their chains
and models as well as the MSA. The returned values are a NamedTuple with the fields
`pdb_files`, `chains`, `models` and `msa`.
"""
function create_alpha_fold_inputs(cluster_folder, ref_pdb, ref_chain, ref_model, 
        pdb_files, chains, models)
    paths = create_pdb_lists(ref_pdb, ref_chain, ref_model, pdb_files, chains, models)
    # create the input MSA
    save_sequences(cluster_folder, paths.pdb_files, chains=paths.chains, models=paths.models)
    msa = align_sequences(joinpath(cluster_folder, "sequences.fasta"))
    cleaned_msa = clean_msa(msa)
    write(joinpath(cluster_folder, "sequences.a3m"), cleaned_msa, MIToS.MSA.FASTA)
    # save the template structures
    template_folder = joinpath(cluster_folder, "templates")
    isdir(template_folder) || mkdir(template_folder)
    for (pdb_file, chain, model) in zip(paths.pdb_files[2:end], paths.chains[2:end], paths.models[2:end])
        pdb_name = basename(pdb_file)
        pdb_data = get_residues_and_sequence(pdb_file; chain=chain, model=model)
        MIToS.PDB.write(joinpath(template_folder, pdb_name), pdb_data.residues, MIToS.PDB.PDBFile)
    end
    
    # return the path to the pdb files, and their chains and models as well as the MSA
    (pdb_files=paths.pdb_files, chains=paths.chains, models=paths.models, msa=cleaned_msa)
end
