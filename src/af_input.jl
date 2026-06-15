# These functions create the folder structure defining AlphaFold 2 inputs, including the 
# the MSA, and the templates.
"""
    _read_pdb_chain(file) -> residues or nothing

Read all heavy atoms from a PDB or mmCIF file.
Input : 
- `file` : path to a `.pdb`, `.cif`, or `.mmcif` file.
Output : 
A residue array (MIToS), or `nothing` if the file cannot be read.
"""
function _read_pdb_chain(file::String) # MIToS.PDB.All
    # assume that the chain_code is All if it is not a string
    # occupancyfilter is needed to avoid the duplicated residue warnings with TMalign
    extension=lowercase(last(splitext(basename(file))))
    if extension == ".cif" || extension == ".mmcif"
        # Read the chain specified by chain_code from the mmcif file
        try
            return MIToS.PDB.read_file(
                file,
                MIToS.PDB.MMCIFFile;
                onlyheavy = true,
                occupancyfilter = true,
            )
        catch err
            @error "Error reading the MMCIF file $file: $err"
            return nothing
        end
    else
        try
            return MIToS.PDB.read_file(
                file,
                MIToS.PDB.PDBFile;
                onlyheavy = true,
                occupancyfilter = true,
            )
        catch err
            @error "Error reading the PDB file $file: $err"
            return nothing
        end
    end
end

"""
    _read_pdb_chain(file, chain_code) -> residues or nothing

Read a specific chain from a PDB or mmCIF file.
Input : 
- `file`       : path to a `.pdb`, `.cif`, or `.mmcif` file.
- `chain_code` : chain identifier to extract (e.g. `"A"`).
Output 
A residue array for the requested chain, or `nothing` on error.
"""
function _read_pdb_chain(file::String, chain_code::String)
    res = nothing
    extension=lowercase(last(splitext(basename(file))))
    if extension == ".cif" || extension == ".mmcif"
        try
            # read the whole file
            res = MIToS.PDB.read_file(
                file,
                MIToS.PDB.MMCIFFile;
                onlyheavy = true,
                occupancyfilter = true,
            )
        catch err
            @error "Error reading the MMCIF file $file: $err"
            return nothing
        end

        # Read the chain specified by chain_code from the PDB file
        # occupancyfilter is needed to avoid the duplicated residue warnings with TMalign
    else
        try
            # read the whole file
            res = MIToS.PDB.read_file(
                file,
                MIToS.PDB.PDBFile;
                onlyheavy = true,
                occupancyfilter = true,
            )
            # note that auth chains can be lower case, so test the lowercase one if 
            # the uppercase one is not found. For example, 7ADD has lowercase chains.

        catch err
            @error "Error reading the PDB file $file: $err"
            return nothing
        end
    end

    try
        chains = Set{String}(r.id.chain for r in res)
        if chain_code in chains
            MIToS.PDB.select_residues(
                res;
                model = MIToS.PDB.All,
                chain = chain_code,
                group = MIToS.PDB.All,
                residue = MIToS.PDB.All,
            )
        else
            lowercase_chain_code = lowercase(chain_code)
            if lowercase_chain_code in chains
                return MIToS.PDB.select_residues(
                    res;
                    model = MIToS.PDB.All,
                    chain = lowercase_chain_code,
                    group = MIToS.PDB.All,
                    residue = MIToS.PDB.All,
                )
            end
            @error "The chain $chain_code was not found in the PDB file $file"
            nothing
        end
    catch err
        @error "Error extracting chain $chain_code from the PDB file $file: $err"
        nothing
    end
end

_read_pdb_chain(file::String, ::Type{MIToS.PDB.All}) = _read_pdb_chain(file)

"""
    get_residues_and_sequence(pdb_file; chain=MIToS.PDB.All, model="1")

This function reads a pdb file and returns the residues and the sequence of the chain
specified by the `chain` argument. The `chain` keyword argument can be a string with the 
chain name or the `MIToS.PDB.All` type. However, this function will throw an error if more 
than one chain is selected; so `All` is used only to avoid specifiying a chain when there is 
only one in the file. The `model` argument is used to select the model to read (by default
the first model is selected).
"""
function get_residues_and_sequence(
    pdb_file;
    chain::Union{Type{MIToS.PDB.All},String} = MIToS.PDB.All,
    model::String = "1",
)
    # if something goes wrong, return an empty vector and an empty string
    query_res = MIToS.PDB.PDBResidue[]
    seq = ""

    res = _read_pdb_chain(pdb_file, chain)
    if res !== nothing
        actual_chain = first(res).id.chain # chain can be lowercase
        #Select all the informations inside the pdb 
        query_res = MIToS.PDB.select_residues(
            res;
            model = model,
            chain = actual_chain,
            group = MIToS.PDB.All,
            residue = MIToS.PDB.All,
        )
        if !isempty(query_res)
            sequences = MIToS.PDB.modelled_sequences(query_res)
            @assert !isempty(sequences) "The are no residues in the pdb file: $(pdb_file) (chain: $actual_chain)"
            if isa(chain, String)
                seq = sequences[(model = model, chain = actual_chain)]
            else
                seq_values = collect(values(sequences))
                if length(seq_values) == 1
                    seq = seq_values[1]
                else
                    @warn "The query must be a single chain if the chain argument is `All`: $(pdb_file)"
                end
            end
        else
            @warn "The are no ATOM residues in $(pdb_file) (model: $model, chain: $actual_chain)"
        end
    end

    (residues = query_res, sequence = seq)
end


"""
    save_sequences(cluster_folder::String, pdb_files::Vector{String}; chains, models)

This function uses the `get_residues_and_sequence` function to create a fasta file with 
the sequences of the structures in the cluster. The second arguments should be a vector
with the paths to the pdb files of the structures in the cluster. The `chains` and `models`
arguments are vectors with the chain and model to read for each pdb file. By default, the
first model and all the chains are selected.
"""
function save_sequences(
    cluster_folder::String,
    pdb_files::Vector{String};
    chains = fill(MIToS.PDB.All, length(pdb_files)),
    models::Vector{String} = fill("1", length(pdb_files)),
)
    open(joinpath(cluster_folder, "sequences.fasta"), "w") do file
        for (pdb_file, chain, model) in zip(pdb_files, chains, models)
            pdb_name = basename(pdb_file)
            pdb_data = get_residues_and_sequence(pdb_file; chain = chain, model = model)
            if !isempty(pdb_data.sequence)
                println(file, ">$pdb_name")
                println(file, pdb_data.sequence)
            end
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
        run(pipeline(`$(MAFFT_jll.mafft()) --quiet $sequences`, stdout = aln_path))
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

    (; pdb_files = pdb_abspaths, chains, models)
end

"""
    create_msa_and_templates(cluster_folder, ref_pdb, ref_chain, ref_model, 
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
function create_msa_and_templates(
    cluster_folder,
    ref_pdb,
    ref_chain,
    ref_model,
    pdb_files,
    chains,
    models,
)
    paths = create_pdb_lists(ref_pdb, ref_chain, ref_model, pdb_files, chains, models)
    # create the input MSA
    save_sequences(
        cluster_folder,
        paths.pdb_files,
        chains = paths.chains,
        models = paths.models,
    )
    msa = align_sequences(joinpath(cluster_folder, "sequences.fasta"))
    cleaned_msa = clean_msa(msa)
    write(joinpath(cluster_folder, "sequences.a3m"), cleaned_msa, MIToS.MSA.FASTA)
    # save the template structures
    template_folder = joinpath(cluster_folder, "templates")
    isdir(template_folder) || mkdir(template_folder)
    # AF: PDB files should have only one model.
    # AF: Filenames should be in lowercase.
    for (pdb_file, chain, model) in
        zip(paths.pdb_files[2:end], paths.chains[2:end], paths.models[2:end])
        pdb_name = basename(pdb_file)
        pdb_code, _ = _get_pdb_and_chain(pdb_name)
        pdb_template = joinpath(template_folder, "$(lowercase(pdb_code)).pdb")
        if !isfile(pdb_template)
            # Keep all the chains in the template pdb structures
            pdb_data =
                get_residues_and_sequence(pdb_file; chain = MIToS.PDB.All, model = model)
            # AF: Insertion codes are not supported.
            res = filter!(r -> isnothing(match(r"[^0-9]+", r.id.number)), pdb_data.residues)
            # Save the cleaned pdb file
            MIToS.PDB.write_file(pdb_template, res, MIToS.PDB.PDBFile)
        end
    end

    # return the path to the pdb files, and their chains and models as well as the MSA
    (
        pdb_files = paths.pdb_files,
        chains = paths.chains,
        models = paths.models,
        msa = cleaned_msa,
    )
end
