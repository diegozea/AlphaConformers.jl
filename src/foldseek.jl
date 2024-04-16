# Functions to run Foldseek inside the pipeline

const _M8_COL_NAMES = ["query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits"]

"""
    foldseek_search(pdb_file::AbstractString; db_path::String = get(ENV, "FOLDSEEK_DB_PATH", ""), format_mode::Int = 0)

Searches a given protein structure in the Foldseek database. 

The `pdb_file` is the path to the protein structure file in PDB format. The `db_path` is 
the path to the Foldseek database, which defaults to the `FOLDSEEK_DB_PATH` environment 
variable. If the `FOLDSEEK_DB_PATH` environment variable is not set and `db_path` is not 
given, the function will throw an error. 

The function executes a Foldseek easy-search. It uses `--format-mode` 0 by default, which
returns the results in the standard m8 format. Those results are saved in an output file 
with the same base name as the PDB file but with `_results.m8` appended to it. The function 
returns the path to the Folseek easy-search output file.

If `format_mode` is set to 5, the function will return the path to the folder with the
aligned structures. The aligned structures contains only the Calpha atoms aligned to the
query structure.

The function operates in a temporary directory that is automatically cleaned up afterwards.
"""
function foldseek_search(pdb_file::AbstractString;
    db_path::String=get(ENV, "FOLDSEEK_DB_PATH", ""),
    format_mode::Int=0)
    isempty(db_path) && error("Please set the FOLDSEEK_DB_PATH environment variable or " *
                              "the db_path keyword argument to the path of the Foldseek database.")
    isfile(db_path) || error("The path to the Foldseek database is not a file.")
    mktempdir() do tmp_folder
        # path without extension
        path = first(splitext(abspath(pdb_file)))
        if format_mode == 5
            folder = dirname(path)
            aligned_folder = joinpath(folder, "aligned_structures")
            if !isdir(aligned_folder)
                mkdir(aligned_folder)
            end
            out_file = joinpath(aligned_folder, "aln_")
        else
            out_file = "$(path)_results.m8"
        end
        run(`$(Foldseek_jll.foldseek()) easy-search $pdb_file $db_path $out_file $tmp_folder --format-mode $format_mode`)
        if format_mode == 5
            rm(out_file) # This file is empty when using using Foldseek v8
            return dirname(out_file)
        else
            return out_file
        end
    end
end

# foldseek ------------------------------------------------------------------------------- #

"""
    read_foldseek_search_results(file::AbstractString)

Reads the Foldseek easy-search output file (m8) and returns a DataFrame with the results 
and proper column names.
"""
function read_foldseek_search_results(file::AbstractString)
    DataFrames.DataFrame(CSV.File(file, delim='\t', header=_M8_COL_NAMES))
end

function run_foldseek(pdb_file::AbstractString,
    db_path::String=get(ENV, "FOLDSEEK_DB_PATH", "");
    out_folder::String=dirname(abspath(pdb_file)))
    # get the path to the target database
    isempty(db_path) && error("Please set the FOLDSEEK_DB_PATH environment variable or " *
                              "the db_path keyword argument to the path of the Foldseek database.")

    # if there is more than one database, run the function for each one
    if occursin(',', db_path)
        db_path_vector = String.(split(db_path, ','))
        return run_foldseek(pdb_file, db_path_vector, out_folder=out_folder)
    end

    # if there is only one database, continue
    isfile(db_path) || error("Foldseek database error: $db_path is not a file.")
    db_name = basename(db_path)

    # IO paths
    out_folder_db = joinpath(out_folder, "$(db_name)_results")
    _create_empty_folder(out_folder_db)
    pdb_name = first(splitext(basename(pdb_file))) # filename without extension
    table_file = joinpath(out_folder_db, "$(pdb_name)_results.m8")
    msa_file = joinpath(out_folder_db, "msa.a3m")
    aligned_structures_folder = joinpath(out_folder_db, "aligned_structures")
    _create_empty_folder(aligned_structures_folder)
    cwd = pwd()

    try
        mktempdir() do tmp_folder
            cd(tmp_folder)

            # createdb for the query file
            run(`$(Foldseek_jll.foldseek()) createdb $pdb_file query_db`)

            # run the search using -a to be able to recover the alignment
            # --prefilter-mode 1 to use less RAM when searching the AFDB, needing ~35 Gb
            run(`$(Foldseek_jll.foldseek()) search query_db $db_path results tmp -a --prefilter-mode 1`)

            # convertalis to m8
            run(`$(Foldseek_jll.foldseek()) convertalis query_db $db_path results $table_file`)

            # convertalis to aligned_structures
            prefix = joinpath(aligned_structures_folder, "aln_")
            run(`$(Foldseek_jll.foldseek()) convertalis query_db $db_path results $prefix --format-mode 5`)
            isfile(prefix) && rm(prefix)

            # run result2msa
            run(`$(Foldseek_jll.foldseek()) result2msa query_db $db_path results msa --msa-format-mode 6`)
            # unpack the msa
            run(`$(Foldseek_jll.foldseek()) unpackdb msa msa_output --unpack-suffix a3m --unpack-name-mode 0`)
            if isfile(msa_file)
                @warn "$msa_file already exists. It will be overwritten."
                rm(msa_file)
            end
            isfile("msa_output/0a3m") && cp("msa_output/0a3m", msa_file)
        end
    finally
        cd(cwd)
    end
    if isfile(msa_file)
        # TODO: Fix on MIToS
        # run(pipeline(`sed '/^>/!s/[a-z]//g' $msa_file`, stdout=msa_file))
        cleaned_file = replace(msa_file, ".a3m" => "_cleaned.a3m")
        open(cleaned_file, "w") do out
            open(msa_file, "r") do in
                for line in eachline(in)
                    if startswith(line, '>')
                        println(out, line)
                    else
                        println(out, replace(line, r"[a-z]" => ""))
                    end
                end
            end
        end
        mv(msa_file, msa_file * ".bak"; force=true)
        mv(cleaned_file, msa_file; force=true)
    end

    [(; table_file, msa_file, aligned_structures_folder)]
end

function run_foldseek(pdb_file::AbstractString, db_path::Vector{String};
    out_folder::String=dirname(abspath(pdb_file)))
    map(db_path) do db
        only(run_foldseek(pdb_file, db, out_folder=out_folder))
    end
end

function merge_tables(table_files::Vector{String})
    tables = map(table_files) do file
        table = read_foldseek_search_results(file)
        table.file .= abspath(file)
        table
    end
    merged = DataFrames.vcat(tables...)
    DataFrames.unique!(merged, ["query", "target", "qstart", "qend", "tstart", "tend"])
end

# create_pdb_folder & structural_clustering ---------------------------------------------- #

function _get_paths(table_file)
    table_folder = dirname(table_file)
    (;
        table=table_file,
        msa=joinpath(table_folder, "msa.a3m"),
        structures=joinpath(table_folder, "aligned_structures")
    )
end

# concatenate msas ----------------------------------------------------------------------- #

function _get_seq_and_columns(msa, pos_ref)
    ref = replace(MIToS.MSA.stringsequence(msa, pos_ref), "-" => "")
    col = [col for (col, res) in enumerate(MIToS.MSA.getsequencemapping(msa, pos_ref))
           if res != 0]
    (ref, col)
end

"""
    get_aligned_positions(aln)

Returns the aligned positions in the two sequences of an alignment. The alignment is any
iterable object that returns a tuple with two elements. For example a 
`BioAlignments.PairwiseAlignment` object as the one returned by `BioAlignments.alignment` 
or a vector of tuples where each tuples is a column of the alignment. The gap should be 
represented by the '-' character. This function return a vector of tuples, each tuple 
containing the aligned positions in the two sequences. For example:

```
julia> get_aligned_positions([('A', 'A'), ('C', '-'), ('G', 'G'), ('T', 'T')])
3-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (3, 2)
 (4, 3)

```
"""
function get_aligned_positions(aln)
    pos_a = 0
    pos_b = 0
    positions = Tuple{Int,Int}[]
    for (res_a, res_b) in aln
        if res_a != '-'
            pos_a += 1
        end
        if res_b != '-'
            pos_b += 1
        end
        if res_a != '-' && res_b != '-'
            push!(positions, (pos_a, pos_b))
        end
    end
    positions
end

function _match_columns(msa_a, msa_b, pos_ref_a, pos_ref_b)
    ref_a, col_a = _get_seq_and_columns(msa_a, pos_ref_a)
    ref_b, col_b = _get_seq_and_columns(msa_b, pos_ref_b)
    aln_model = BioAlignments.AffineGapScoreModel(match=6, mismatch=-4,
        gap_open=-2, gap_extend=-1)
    aln = BioAlignments.alignment(BioAlignments.pairalign(
        BioAlignments.GlobalAlignment(), ref_a, ref_b, aln_model))
    positions = get_aligned_positions(aln)
    columns_a = [col_a[pos[1]] for pos in positions]
    columns_b = [col_b[pos[2]] for pos in positions]
    (columns_a, columns_b)
end

function merge_msas(table::DataFrames.DataFrame)
    out_folders = dirname.(unique(table.file))
    msas = [read(joinpath(folder, "msa.a3m"), MIToS.MSA.FASTA, generatemapping=true)
            for folder in out_folders]
    # Return the MSA if there is only one
    length(msas) == 1 && return only(msas)
    # Otherwise, concatenate the multiple MSAs using the query sequence as reference
    msa_a = msas[1]
    for msa_b in msas[2:end]
        cols_a, cols_b = _match_columns(msa_a, msa_b, 1, 1)
        msa_a = MIToS.MSA.join(msa_a, msa_b[2:end, :], cols_a, cols_b; axis=2)
    end
    # Clean the sequence names by deleting the MSA number at the beginning that 
    # was added by the join function and the sequence data Foldseek adds at the end.
    cleaned_names = String[replace(first(split(name)), r"^[1-9]+_" => "")
                           for name in MIToS.MSA.sequencename_iterator(msa_a)]
    MIToS.MSA.rename_sequences!(msa_a, cleaned_names)
end

# structure_distances -------------------------------------------------------------------- #

function _aligned_pdb_file(query, target, aligned_structures_folder)
    file_name = "aln_$(query)_$(target).pdb"
    joinpath(aligned_structures_folder, file_name)
end

function _get_aligned_residues(query, target, paths, tstart, tend)
    # msa = read(paths.msa, MIToS.MSA.FASTA)
    pdb_file = _aligned_pdb_file(query, target, paths)
    res = read(pdb_file, MIToS.PDB.PDBFile)
    # NOTE: Foldseek creates a single chain for the aligned structures: A
    # NOTE: tstart and tend are the residue numbers in the Foldseek-aligned structure
    res[tstart:tend]
end

# NOTE: This function only works before adding external structures into the
# Foldseek results, as it relies on the structures generated by Foldseek.
function _get_aligned_structures(foldseek_results::DataFrames.DataFrame)
    structures = OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}}()
    for row in eachrow(foldseek_results)
        paths = _get_paths(row.file)
        try
            structures[row.target] = _get_aligned_residues(row.query, row.target, 
                paths.structures, row.tstart, row.tend)
        catch err
            @error "Error ($err) getting structures from row: $row"
            continue
        end
    end
    structures
end

# add known conformations ---------------------------------------------------------------- #
#
# This functions add strcutural information from known conformations added to the Foldseek
# results after running `add_known_conformations!`.
#
# First, I need a function to perform the structural alignment between the conformations
# using MIToS.

# To test I can use:
# conformation_a = structures["4F4J.pdb_A"]; conformation_b = structures["1EX7.pdb"]

"""
    coverage_and_identity(aln)

It returns the alignment coverage and the percentage of identity relative to the 
first sequence. The alignment is any iterable object that returns a tuple with two 
elements. For example a `BioAlignments.PairwiseAlignment` object as the one returned 
by `BioAlignments.alignment` or a vector of tuples where each tuples is a column of 
the alignment. The gap should be represented by the '-' character.

```julia
julia> coverage_and_identity(collect(zip("AA-AAA", "A-CCAA")))
(0.8, 0.75)

```

i.e., 4 out of 5 residues match in the first sequence, and among these, 3 out of 4 
residues are identical
"""
function coverage_and_identity(aln)
    total_a = 0
    matched_a = 0
    identities = 0
    for (res_a, res_b) in aln
        if res_a != '-'
            total_a += 1
            if res_b != '-'
                matched_a += 1
                if res_a == res_b
                    identities += 1
                end
            end
        end
    end
    (matched_a / total_a, identities / matched_a)
end

function structural_alignment(conformation_a, conformation_b,
    aln_type=BioAlignments.OverlapAlignment(),
    aln_model=BioAlignments.AffineGapScoreModel(match=6, mismatch=-4, gap_open=-2,
        gap_extend=-1))
    # only keep residues with 'CA' atom
    clean_a = filter(res -> !isempty(MIToS.PDB.findatoms(res, "CA")), conformation_a)
    clean_b = filter(res -> !isempty(MIToS.PDB.findatoms(res, "CA")), conformation_b)
    len_a = length(clean_a)
    len_b = length(clean_b)
    if len_a == 0 || len_b == 0
        @warn "One of the structures has no 'CA' atoms: $len_a, $len_b"
        return nothing
    end
    try
       # get the sequences
        seq_a = first(values(MIToS.PDB.modelled_sequences(clean_a)))
        seq_b = first(values(MIToS.PDB.modelled_sequences(clean_b)))
        # align the sequences
        aln = BioAlignments.alignment(BioAlignments.pairalign(aln_type, seq_a, seq_b, aln_model))
        # get the stats
        coverage, identity = coverage_and_identity(aln)
        # get the aligned residues
        matches = get_aligned_positions(aln)
        # structural superposition of the aligned residues
        aligned_a, aligned_b, rmsd = MIToS.PDB.superimpose(clean_a, clean_b, matches)
        return (aligned_a, aligned_b, matches, rmsd, coverage, identity)
    catch err
        @error "Error in the structural alignment: $err"
        return nothing
    end
end

function _read_pdb(pdb_file, chain)
    if isfile(pdb_file)
        MIToS.PDB.read(pdb_file, MIToS.PDB.PDBFile,
            chain=chain, model="1", onlyheavy=true, occupancyfilter=true)
    else
        @warn "The file $pdb_file does not exist."
        nothing
    end
end

# adds the best match to the structures dictionary
function find_best_match!(
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}},
    rmsds::Dict{Tuple{String,String},Float64},
    target::String,
    target2uniprot::Dict{String,String},
    uniprot2targets::Dict{String,Vector{String}};
    pdb_folder::Union{String,Nothing}=nothing,
    min_coverage::Float64=0.75,
    min_identity::Float64=0.95)
    # get the PDB code and chain from the target name
    pdb, chain = AlphaConformers._get_pdb_and_chain(target)
    # read the PDB file
    conformation_b = if pdb_folder === nothing
        mktempdir() do tmp_folder
            @info "Downloading $pdb"
            pdb_file = joinpath(tmp_folder, "$pdb.pdb.gz")
            MIToS.PDB.downloadpdb(pdb; format=MIToS.PDB.PDBFile, filename=pdb_file)
            _read_pdb(pdb_file, chain)
        end
    else
        _read_pdb(joinpath(pdb_folder, uppercase(pdb) * ".pdb"), chain)
    end
    if conformation_b === nothing
        return nothing
    end
    # get the structures that comes from the same UniProt entry
    targets = uniprot2targets[target2uniprot[target]]
    # perform all the structural alignments
    alignment_list = []
    map(targets) do target_a
        conformation_a = structures[target_a]
        sorted_targets = sort([target, target_a])
        key = (sorted_targets[1], sorted_targets[2])
        aln = structural_alignment(conformation_a, conformation_b)
        if aln !== nothing
            push!(alignment_list, key => aln)
        end
    end
    if isempty(alignment_list)
        return nothing
    end
    alignments = OrderedCollections.OrderedDict(alignment_list)
    # get the best match
    df = DataFrames.DataFrame(values(alignments), ["aligned_a", "aligned_b", "matches",
        "rmsd", "coverage", "identity"])
    df[!, "keys"] .= keys(alignments)
    filter!(row -> row.coverage > min_coverage && row.identity > min_identity, df)
    if !isempty(df)
        # if there is a match, store the rmsd values
        for (key, aln) in alignments
            rmsds[key] = aln[4]
        end
        # select the best match
        sort!(df, ["identity", "coverage", "rmsd"], rev=[true, true, false])
        best_alignment = first(df) # return a DataFrameRow containing the best match
        # store the best match in the structures dictionary
        aligned_positions = [m[2] for m in best_alignment.matches]
        structures[target] = best_alignment.aligned_b[aligned_positions]
        # return the best match
        best_alignment
    else
        nothing
    end
end

function process_known_conformations!(
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}},
    expanded_table::DataFrames.DataFrame,
    target2uniprot::Dict{String,String},
    uniprot2targets::Dict{String,Vector{String}};
    pdb_folder::Union{String,Nothing}=nothing,)
    rmsds = Dict{Tuple{String,String},Float64}()
    for row in eachrow(expanded_table)
        if ismissing(row.query)
            target = row.target
            best_match = find_best_match!(structures, rmsds, target,
                target2uniprot, uniprot2targets, pdb_folder=pdb_folder)
            if best_match !== nothing
                row.query = only([q for q in best_match.keys if q != target])
            end
        end
    end
    filter!(row -> !ismissing(row.query), expanded_table)
    rmsds
end

function fill_rmsds!(rmsds::Dict{Tuple{String,String},Float64},
    table::DataFrames.DataFrame,
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}})
    n = DataFrames.nrow(table)
    ij_pairs = collect(Combinatorics.combinations(1:n, 2))
    ProgressMeter.@showprogress Distributed.@distributed for (i, j) in ij_pairs
        row_i = table[i, :]
        row_j = table[j, :]
        sorted_ids = sort([row_i.target, row_j.target])
        key = (sorted_ids[1], sorted_ids[2])
        if haskey(rmsds, key)
            continue
        end
        struct_a = structures[row_i.target]
        struct_b = structures[row_j.target]
        aln = structural_alignment(struct_a, struct_b,
            BioAlignments.GlobalAlignment(),
            BioAlignments.AffineGapScoreModel(BioAlignments.BLOSUM62,
                gap_open=-10, gap_extend=-1))
        if aln !== nothing
            rmsds[key] = aln[4]
        else
            rmsds[key] = NaN
        end
    end
    rmsds
end

function get_rmsd_matrix(rmsds::Dict{Tuple{String,String},Float64}, targets::Set{String})
    target_a = String[]
    target_b = String[]
    rmsd = Float64[]
    for (key, value) in rmsds
        key_a = key[1]
        key_b = key[2]
        if key_a in targets && key_b in targets
            push!(target_a, key_a)
            push!(target_b, key_b)
            push!(rmsd, value)
        end
    end
    df = DataFrames.DataFrame(; target_a, target_b, rmsd)
    sort!(df, [:target_a, :target_b])
    mat = PairwiseListMatrices.from_table(df, false, diagonalvalue=0.0)
    # delete rows and columns with NaN values
    nan_col = vec(any(isnan, mat, dims=1))
    nan_row = vec(any(isnan, mat, dims=2))
    mat[.!(nan_row), .!(nan_col)]
end

function cluster_structures(rmsds::Dict{Tuple{String,String},Float64}, targets::Set{String})
    rmsd_mat = get_rmsd_matrix(rmsds, targets)
    unsorted_names = names(rmsd_mat, 1)
    clustering_result = Clustering.hclust(rmsd_mat, linkage=:complete, branchorder=:optimal)
    structure_names = unsorted_names[clustering_result.order]
    structure_names, clustering_result
end

function get_target2sequence(expanded_table, msa)
    seqnames = Set{String}(MIToS.MSA.sequencename_iterator(msa))
    target2sequence = Dict{String,String}()
    for row in eachrow(expanded_table)
        seqname = ismissing(row.evalue) ? row.query : row.target
        if seqname in seqnames
            target2sequence[row.target] = seqname
        else
            @warn "The sequence $seqname is not in the MSA."
        end
    end
    target2sequence
end

function get_cluster2targets(targets, clusters)
    cluster2targets = OrderedCollections.OrderedDict{Int,Vector{String}}()
    cluster_numbers = unique(clusters)
    for cluster in cluster_numbers
        cluster2targets[cluster] = targets[clusters.==cluster]
    end
    cluster2targets
end

function get_cluster2seqnames(cluster2targets, target2sequence)
    cluster2seqnames = OrderedCollections.OrderedDict{Int,Vector{String}}()
    for (cluster, targets) in cluster2targets
        cluster2seqnames[cluster] = unique(target2sequence[target] for target in targets
                                           if haskey(target2sequence, target))
    end
    cluster2seqnames
end

function get_cluster2msa(msa, cluster2seqnames)
    query_name = first(MIToS.MSA.sequencename_iterator(msa))
    cluster2msa = OrderedCollections.OrderedDict{Int,MIToS.MSA.AnnotatedMultipleSequenceAlignment}()
    for (cluster, seqnames) in cluster2seqnames
        seqnames = [query_name; seqnames]
        cluster2msa[cluster] = msa[seqnames, :]
    end
    cluster2msa
end

function get_cluster2structures(structures, cluster2targets)
    cluster2structures = OrderedCollections.OrderedDict{Int,Dict{String,Vector{MIToS.PDB.PDBResidue}}}()
    for (cluster, targets) in cluster2targets
        cluster2structures[cluster] = Dict{String,Vector{MIToS.PDB.PDBResidue}}(
            target => structures[target] for target in targets)
    end
    cluster2structures
end

function create_template_clusters(rmsds::Dict{Tuple{String,String},Float64},
    expanded_table::DataFrames.DataFrame,
    msa::MIToS.MSA.AnnotatedMultipleSequenceAlignment,
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}})
    target2sequence = AlphaConformers.get_target2sequence(expanded_table, msa)
    targets = Set{String}(expanded_table.target)
    targets, clustering_result = cluster_structures(rmsds, targets)
    large_clusters = Clustering.cutree(clustering_result, h=1.0)
    small_clusters = Clustering.cutree(clustering_result, h=0.5)
    large_cluster2targets = AlphaConformers.get_cluster2targets(targets, large_clusters)
    small_cluster2targets = AlphaConformers.get_cluster2targets(targets, small_clusters)
    large_cl2seq = AlphaConformers.get_cluster2seqnames(large_cluster2targets, target2sequence)
    small_cl2seq = AlphaConformers.get_cluster2seqnames(small_cluster2targets, target2sequence)
    large_cl2msa = AlphaConformers.get_cluster2msa(msa, large_cl2seq)
    small_cl2pdb = AlphaConformers.get_cluster2structures(structures, small_cluster2targets)
    large_small_pairs = zip(large_clusters, small_clusters) |> unique |> sort
    (large_small_pairs, large_cl2msa, small_cl2pdb)
end

function create_folder_structure(large_small_pairs::Vector{Tuple{Int,Int}},
    large_cl2msa::OrderedCollections.OrderedDict{Int,MIToS.MSA.AnnotatedMultipleSequenceAlignment},
    small_cl2pdb::OrderedCollections.OrderedDict{Int,Dict{String,Vector{MIToS.PDB.PDBResidue}}};
    out_folder::String=mktempdir())
    for (large, small) in large_small_pairs
        # MSA
        cluster_folder = mkdir(joinpath(out_folder, "cluster_$(large)_$(small)"))
        msa_file = joinpath(cluster_folder, "sequences.a3m")
        write(msa_file, large_cl2msa[large], MIToS.MSA.FASTA)
        # Structures
        cluster_template_folder = mkdir(joinpath(cluster_folder, "templates"))
        for (target, structure) in small_cl2pdb[small]
            pdb_file = joinpath(cluster_template_folder,
                replace(target, ".pdb" => "") * ".pdb")
            write(pdb_file, structure, MIToS.PDB.PDBFile)
        end
    end
    out_folder
end

function alphaconformers(input_pdb, pdb_db, alphafold_db, pdb_folder, out_folder; 
        evalue_cutoff::Float64=1e-5)
    
    @info "Running Foldseek"
    output = run_foldseek(input_pdb, "$pdb_db,$alphafold_db", out_folder=out_folder)
    merged_table = merge_tables([output[1].table_file, output[2].table_file])
    filter!(row -> row.evalue < evalue_cutoff, merged_table)
    merged_msa = merge_msas(merged_table)
    
    @info "Getting the aligned structures"
    structures = _get_aligned_structures(merged_table)
    
    @info "Adding known conformations"
    sifts_uniprot_mapping = get_uniprot_mapping()
    target2uniprot, expanded_table = add_known_conformations!(deepcopy(merged_table), sifts_uniprot_mapping)
    uniprot2targets = get_uniprot2targets(target2uniprot, expanded_table)
    
    @info "Measuring RMSDs"
    rmsds = process_known_conformations!(structures, expanded_table, target2uniprot, uniprot2targets, pdb_folder=pdb_folder)
    fill_rmsds!(rmsds, expanded_table, structures)
    
    @info "Clustering structures"
    large_small_pairs, large_cl2msa, small_cl2pdb = create_template_clusters(rmsds, expanded_table, merged_msa, structures)
    create_folder_structure(large_small_pairs, large_cl2msa, small_cl2pdb, out_folder=out_folder)
end


# USAGE EXAMPLE
# =============
#=


#using Revise, AlphaConformers, DataFrames, MIToS; 
#input_pdb = joinpath("/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/AlphaConformers/test", "data", "1EX6_B.pdb"); 
#test_db = joinpath("/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/AlphaConformers/test", "data", "test_db", "test_db"); 
#test_db_2 = joinpath("/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/AlphaConformers/test", "data", "test_db_2", "test_db_2"); 
#temp_dir = mktempdir()

using Revise, AlphaConformers, DataFrames, MIToS; input_pdb = joinpath("/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/AlphaConformers/test", "data", "1EX6_B.pdb"); test_db = joinpath("/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/AlphaConformers/test", "data", "test_db", "test_db"); test_db_2 = joinpath("/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/AlphaConformers/test", "data", "test_db_2", "test_db_2"); temp_dir = mktempdir()

AlphaConformers.alphaconformers(input_pdb, test_db, test_db_2, temp_dir)

AlphaConformers.run_alphafold(temp_dir, colabfold_path="/opt/alphafold/runcolabfold.py")



output = run_foldseek(input_pdb, "$test_db,$test_db_2", out_folder=temp_dir)
merged_table = merge_tables([output[1].table_file, output[2].table_file])
merged_msa = merge_msas(merged_table)
structures = AlphaConformers._get_aligned_structures(merged_table)

sifts_uniprot_mapping = get_uniprot_mapping()

target2uniprot, expanded_table = AlphaConformers.add_known_conformations!(deepcopy(merged_table), sifts_uniprot_mapping)

uniprot2targets = AlphaConformers.get_uniprot2targets(target2uniprot, expanded_table)

rmsds = AlphaConformers.process_known_conformations!(structures, expanded_table, target2uniprot, uniprot2targets)

AlphaConformers.fill_rmsds!(rmsds, expanded_table, structures)

structure_names, clustering_result = AlphaConformers.cluster_structures(rmsds)

clusters = Clustering.cutree(clustering_result, h=1.0)

large_small_pairs, large_cl2msa, small_cl2pdb = AlphaConformers.create_template_clusters(rmsds, expanded_table, merged_msa, structures)

out_folder = AlphaConformers.create_folder_structure(large_small_pairs, large_cl2msa, small_cl2pdb)

# [ ] TODO : Filter Foldseek results based on E-values.

# [ ] TODO: I should make `add_known_conformations!` to return a directory from query/target to 
# UniProt (or the other way around). Then, I can also create the other one from that one.
# The idea is to be able to align the known conformations to the Foldseek results using 
# sequence alignment. I will need to match the PDB sequences to the MSA sequences, and add 
# the aligned residues into the structures.

# [ ] WARNING: What does it happen if I have a single domain, and one of the target has two 
# identical domains? Since structures use the target name as a dictionary key, I think I 
# will overwritte the first domain with the second one. I should think what to do in that case.

# list_known_conformations(merged_table[4:4, :], sifts_uniprot_mapping) # It works with AFDB
# list_known_conformations(merged_table, sifts_uniprot_mapping)
# AlphaConformers.add_known_conformations!(deepcopy(merged_table), sifts_uniprot_mapping)
# AlphaConformers._get_pdb_and_chain("1EX6_B.pdb")
# delete_query_from_target!(deepcopy(merged_table), sifts_uniprot_mapping, "1EX6", "B")
# # all nice until here, now it comes the things we need to change in af_input.jl :
# # create_pdb_folder
# # structural_clustering

=#


# TODO: create_alpha_fold_inputs should change to use the new Foldseek results
# The next step is to add: 
# [ ] 1. Add the external structures to the Foldseek results, for this, we should align the 
# conformers to the know aligned conformation; allowing for adding possible fragmented 
# structures into the MSA and Foldseek tables.
#    # This is in the TODO before this one. Finally, I have decided I can use sequence 
#    # alignment to add the known conformations to the Foldseek results. But, I should 
#    # keep also track of the aligned structural residues for later use.
# [ ] 3. Use the MSA matching to calculate the RMSD and clusterize.
# [ ] 4. Filter base on RMSD (<20) and coverage.
# [ ] 2. Once we have all the structures, we can create the input files for AlphaFold.

#=
function structure_distances(foldseek_results::DataFrames.DataFrame)
    targets = foldseek_results.target
    n_targets = length(targets) # also nrows
    distances = zeros(Float16, n_targets, n_targets)
    for i in 1:(n_targets-1)
        row_i = foldseek_results[i, :]
        paths_i = _get_paths(row_i.file)
        pdb_i = _get_aligned_residues(row_i.query, row_i.target, paths_i.structures, 
            row_i.tstart, row_i.tend)
        for j in (i+1):n_targets
            row_j = foldseek_results[j, :]
            paths_j = _get_paths(row_j.file)
            pdb_j = _get_aligned_residues(row_j.query, row_j.target, paths_j.structures, 
                row_j.tstart, row_j.tend)
            @show MIToS.PDB.rmsd(pdb_i, pdb_j)
            distances[i, j] = distances[j, i] = MIToS.PDB.rmsd(pdb_i, pdb_j) 
        end
    end
end
=#