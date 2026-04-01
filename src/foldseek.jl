# Functions to run Foldseek inside the pipeline

const _M8_COL_NAMES = ["query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits","qtmscore", "ttmscore", "alntmscore", "rmsd", "prob"]




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
    format_output::String="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore"
)
    isempty(db_path) && error("Please set the FOLDSEEK_DB_PATH environment variable or " *
                              "the db_path keyword argument to the path of the Foldseek database.")
    isfile(db_path) || error("The path to the Foldseek database is not a file.")
    mktempdir() do tmp_folder
        # path without extension
        path = first(splitext(abspath(pdb_file)))
        #=
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
        =#
        out_file = "$(path)_results.m8"
        run(`$(Foldseek_jll.foldseek()) easy-search $pdb_file $db_path $out_file $tmp_folder --format-output $format_output`)
        #=
        if format_mode == 5
            rm(out_file) # This file is empty when using using Foldseek v8
            return dirname(out_file)
        else
            return out_file
        end
        =#
        return out_file
    end
end

# foldseek ------------------------------------------------------------------------------- #

"""
    read_foldseek_search_results(file::AbstractString)

Reads the Foldseek easy-search output file (m8) and returns a DataFrame with the results 
and proper column names.
"""
function read_foldseek_search_results(file::AbstractString; colonnes::Vector{String}=_M8_COL_NAMES)
    """
    lines = readlines(file)
    # Garder seulement les lignes qui ne commencent pas par "@SQ"
    filtered_lines = filter(line -> !startswith(line, "@SQ"), lines)
    # Écrire dans un fichier temporaire
    temp_file = tempname()
    open(temp_file, "w") do io
        for line in filtered_lines
            println(io, line)
        end
    end
    """
    # Charger dans un DataFrame
    df = DataFrames.DataFrame(CSV.File(file, delim='\t', header=colonnes))
    return df
    #rm(temp_file; force=true)
end

function run_foldseek(pdb_file::AbstractString,
    n_threads::Int,
    db_path::String=get(ENV, "FOLDSEEK_DB_PATH", "");
    out_folder::String=dirname(abspath(pdb_file)),
    filtrage::Bool=true,
    format_output::String="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qtmscore,ttmscore,alntmscore,rmsd,prob")
    # get the path to the target database
    isempty(db_path) && error("Please set the FOLDSEEK_DB_PATH environment variable or " *
                              "the db_path keyword argument to the path of the Foldseek database.")

    # if there is more than one database, run the function for each one
    if occursin(',', db_path)
        db_path_vector = String.(split(db_path, ','))
        return run_foldseek(pdb_file, n_threads, db_path_vector, out_folder=out_folder)
    end

    # if there is only one database, continue
    isfile(db_path) || error("Foldseek database error: $db_path is not a file.")
    db_name = basename(db_path)
    @show db_name
    # IO paths
    out_folder_db = joinpath(out_folder, "$(db_name)_results")
    @show out_folder_db
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
            run(`$(Foldseek_jll.foldseek()) version`)

            # createdb for the query file
            run(pipeline(`$(Foldseek_jll.foldseek()) createdb $pdb_file query_db --threads $n_threads`))

            # run the search using -a to be able to recover the alignment
            # --prefilter-mode 1 to use less RAM when searching the AFDB, needing ~35 Gb
            if filtrage 
                run(pipeline(`$(Foldseek_jll.foldseek()) search query_db $db_path results tmp -a -s 10 --max-seqs 1000 -e 10 --prefilter-mode 1 --threads $n_threads`, stdout=joinpath(out_folder_db,"output"), stderr=joinpath(out_folder_db,"error")))
            else 
                run(pipeline(`$(Foldseek_jll.foldseek()) search query_db $db_path results tmp -a -s 1 --max-seqs 1000000 -e inf --prefilter-mode 0 --threads $n_threads`, stdout=joinpath(out_folder_db,"output"), stderr=joinpath(out_folder_db,"error")))
            end
            # convertalis to m8
            run(pipeline(`$(Foldseek_jll.foldseek()) convertalis query_db $db_path results $table_file --format-output $format_output --exact-tmscore 1 --threads $n_threads`))

            # convertalis to aligned_structures
            prefix = joinpath(aligned_structures_folder, "aln_")
            run(pipeline(`$(Foldseek_jll.foldseek()) convertalis query_db $db_path results $prefix --format-mode 5 --threads $n_threads`))
            isfile(prefix) && rm(prefix, recursive=true)

            # run result2msa
            run(pipeline(`$(Foldseek_jll.foldseek()) result2msa query_db $db_path results msa --msa-format-mode 6 --threads $n_threads`))
            # unpack the msa
            run(pipeline(`$(Foldseek_jll.foldseek()) unpackdb msa msa_output --unpack-suffix a3m --unpack-name-mode 0 --threads $n_threads`))
            if isfile(msa_file)
                @warn "$msa_file already exists. It will be overwritten."
                rm(msa_file, recursive=true)
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

function run_foldseek(pdb_file::AbstractString,n_threads::Int,db_path::Vector{String};
    out_folder::String=dirname(abspath(pdb_file)))
    map(db_path) do db
        only(run_foldseek(pdb_file, n_threads, db, out_folder=out_folder))
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

function _seq_name_to_key(sequence_name)
    fields = first(split(sequence_name, "\t"), 3)
    (fields[1], parse(Int, fields[2]), parse(Float64, fields[3]))
end

function _find_duplicates(lst)
    seen = Set()
    duplicates = Set()
    for x in lst
        if x in seen
            push!(duplicates, x)
        else
            push!(seen, x)
        end
    end
    println(length(duplicates))
    println(duplicates)
    return duplicates
end

function merge_msas(table)
    if table isa DataFrames.DataFrame
        println("table ",size(table))
        out_folders = dirname.(unique(table.file))
        starts = Set((row.target, row.bits, row.fident) for row in eachrow(table))
        targets = Set(row.target for row in eachrow(table))
        @show "starts ",first(starts,20)
        @show "targets ",first(targets,5)
        
    else
        out_folders = table
    end
    
    msas = [MIToS.MSA.read_file(joinpath(folder, "msa.a3m"), MIToS.MSA.FASTA, generatemapping=true)
            for folder in out_folders]
    # Select only the matched sequences by using bits and fident to identify the matches.
    # Apply this filter only when there are duplicated names to prevent losing elements 
    # due to numerical differences in comparisons.
    @show "shape msas ",size(msas)
    if size(msas)[1] == 0
        @error "No MSAs found in the provided table."
        return nothing
    end
    for i in eachindex(msas)
        msa = msas[i]
        @show " msa ",size(msa)
        seqnames = MIToS.MSA.sequencenames(msa)[2:end]
        @show "shape seqnames ",length(seqnames)
        if abs(length(seqnames) - size(msa, 1)) != 1
            @info "Avertissement : le nombre de noms de séquences ne correspond pas au nombre de lignes de msa !"
            while length(seqnames) < (size(msa, 1)-1)
                push!(seqnames, "Unnamed_" * string(length(seqnames) + 1))
            end
        end
        @show "shape selected ",length(seqnames)
        @show "Nombre unique d'IDs :", length(unique(seqnames)) 
        @show "shape msa ", size(msa)
        msa_targets = [first(split(seqname, "\t")) for seqname in seqnames]
        @show "msa_target ",length(msa_targets)
        @show length(unique(msa_targets))
        duplicated_msa_targets = _find_duplicates(msa_targets)
        @show duplicated_msa_targets
        # Create a selection vector, starting with `true` to keep the first sequence
        selected = Bool[]
        selected = trues(length(seqnames)+1)
        names = Dict{Tuple{String, Int, Float64}, Int}()  # name → (key, index)
        index=1
        for seqname in seqnames
            index=index+1
            name = first(split(seqname, "\t"))
            if occursin("Unnamed", name)
                 selected[index] = false
            elseif name in duplicated_msa_targets
                key = _seq_name_to_key(seqname)  # ("1QGP.pdb_A", 56, 0.203)
                

                if key in starts
                    # Vérifie si une entrée du même nom existe déjà dans la liste `names`
                    found=nothing
                    for kv in names
                        
                        if kv[1][1] == key[1]
                            
                            found = kv
                        end
                    end

                    #found = findfirst(kv -> kv[1][1] == key[1], names)
                    if found !== nothing
                        existing_key = found[1]
                        existing_index = found[2]
                        
                        # Comparaison sur la valeur du score (3e élément)
                        if key[3] < existing_key[3]
                            
                            # Trouver la position de l'entrée à désélectionner
                            selected[existing_index] = false
                            selected[index] = true
                            names[key] = index
                            delete!(names, existing_key)  # key est la clé, pas la valeur

                        else
                            selected[index] = false
                            
                        end
                    else
                        # Première fois qu'on voit ce nom
                        names[key] = index
                        selected[index] = true
                        
                    end
                else
                    selected[index] = false
                    
                end
            else
                selected[index] = name in targets

            end
        end  
        @show "shape selected ",length(selected)  
        @show "shape msa ", size(msa)
        @show "Nombre de true dans selected : ", count(selected)  # Vérifie combien de séquences sont sélectionnéess
        msas[i] = msa[selected, :]
        @show first(msa,5)
    end
    @show "shape msas ",length(msas)
    # Return the MSA if there is only one
    if length(msas) == 1
        @info "1 msa"
        return only(msas) 
    end
    #length(msas) == 1 && return only(msas)
    # Otherwise, concatenate the multiple MSAs using the query sequence as reference
    msa_a = msas[1]
    for msa_b in msas[2:end]
        cols_a, cols_b = _match_columns(msa_a, msa_b, 1, 1)
        msa_a = MIToS.MSA.join_msas(msa_a, msa_b[2:end, :], cols_a, cols_b; axis=2)
    end
    @show "msa_a new ",first(msa_a,2)
    # Clean the sequence names by deleting the MSA number at the beginning that 
    # was added by the join function and the sequence data Foldseek adds at the end.
    msa_a_seq=collect(MIToS.MSA.sequencename_iterator(msa_a))
    @show first(msa_a_seq,5)
    cleaned_names=String[]
    for name in msa_a_seq
        all_name=split(name,"\t")[1]
        push!(cleaned_names, join(split(all_name, "_")[2:end], "_"))
    end 
    @show first(cleaned_names,5)
    @show length(cleaned_names)
    #cleaned_names = String[replace(first(split(name)), r"^[1-9]+_" => "")
    #                  for name in MIToS.MSA.sequencename_iterator(msa_a)]
    @show length(unique(cleaned_names))
    if length(cleaned_names) != length(unique(cleaned_names))
        @info "Noms dupliqués détectés → suppression des doublons"

        seen = Set{String}()
        selected = trues(length(cleaned_names))

        for (i, name) in enumerate(cleaned_names)
            if name in seen
                selected[i] = false   # on supprime le doublon
            else
                push!(seen, name)
            end
        end

        # Filtrer le MSA
        msa_a = msa_a[selected, :]

        # Mettre à jour les noms nettoyés
        cleaned_names = cleaned_names[selected]
    end
  
    @assert size(msa_a, 1) == length(cleaned_names)

    # Renommer les séquences avec les nouveaux noms uniques
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
    res = MIToS.PDB.read_file(pdb_file, MIToS.PDB.PDBFile)
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
        seqs_a = MIToS.PDB.modelled_sequences(clean_a)
        seqs_b = MIToS.PDB.modelled_sequences(clean_b)

        if isempty(seqs_a) || isempty(seqs_b)
            @warn "No modelled sequences could be extracted (empty sequences)"
            return nothing
        end

        seq_a = first(values(seqs_a))
        seq_b = first(values(seqs_b))
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
        MIToS.PDB.read_file(pdb_file, MIToS.PDB.PDBFile,
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
    conformation_b = if pdb_folder === nothing #link with pdb databases 
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
        if conformation_a === nothing
            @warn "Structure for target $target_a is missing."
            return
        end
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


    names = unique(vcat(df.aligned_a, df.aligned_b))
    N = length(names)

    for name in names
        count_a = count(==(name), df.aligned_a)
        count_b = count(==(name), df.aligned_b)
        total = count_a + count_b
        if total != N - 1
            @warn "Le nom $name apparaît $total fois (attendu : $(N-1))"
        end
    end
    #filter!(row -> !ismissing(row.coverage) > min_coverage && !ismissing(row.identity) > min_identity, df)
    #filter!(row ->  !ismissing(row.identity) > min_identity, df)
    
    if !isempty(df)
        # if there is a match, store the rmsd values
        for (key, aln) in alignments
            if aln[4] === nothing
                return nothing
            else 
                rmsds[key] = aln[4]
            end
        end
        # select the best match
        #sort!(df, ["identity", "coverage", "rmsd"], rev=[true, true, false])
        sort!(df, ["rmsd"])
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
    @show first(expanded_table,10)
    for row in eachrow(expanded_table)
        if ismissing(row.query)
            target = row.target
            try 
                best_match = find_best_match!(structures, rmsds, target,
                    target2uniprot, uniprot2targets, pdb_folder=pdb_folder)
                if best_match !== nothing
                    row.query = only([q for q in best_match.keys if q != target])
                end
            catch e
                @error "Error processing target $target"
                continue
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
    ij_pairs = Tuple{Int, Int}[]
    for i in 1:(n-1)
        for j in i+1:n
            push!(ij_pairs, (i, j))
        end
    end
    # ProgressMeter.@showprogress Distributed.@distributed for (i, j) in ij_pairs
    for (i, j) in ij_pairs
        row_i = table[i, :]
        row_j = table[j, :]
        sorted_ids = sort([row_i.target, row_j.target])
        key = (sorted_ids[1], sorted_ids[2])
        if haskey(rmsds, key)
            continue
        end
        struct_a = structures[row_i.target]
        struct_b = structures[row_j.target]
        if struct_a === nothing || struct_b === nothing
            @warn "Structure for target $(row_i.target) or $(row_j.target) is missing."
            rmsds[key] = nothing
            continue
        end
        aln = structural_alignment(struct_a, struct_b,
            BioAlignments.GlobalAlignment(),
            BioAlignments.AffineGapScoreModel(BioAlignments.BLOSUM62,
                gap_open=-10, gap_extend=-1))
        if aln !== nothing
            rmsds[key] = aln[4]
        else
            rmsds[key] = nothing
        end
    end
    rmsds
end

function get_rmsd_matrix(rmsds::Dict{Tuple{String,String},Float64}, targets::Set{String})
    target_a = String[]
    target_b = String[]
    rmsd = Float64[]
    @show size(targets)
    for (key, value) in rmsds
        key_a = key[1]
        key_b = key[2]
        if key_a in targets && key_b in targets
            push!(target_a, key_a)
            push!(target_b, key_b)
            if value===nothing
                push!(rmsd, NaN)
            else 
                push!(rmsd, value)
            end
        end
    end
    df = DataFrames.DataFrame(; target_a, target_b, rmsd)
    @show size(df)
    sort!(df, [:target_a, :target_b])
    for col in names(df)
        col_data = df[!, col]

        n_missing = count(ismissing, col_data)
        n_nothing = count(x -> x === nothing, col_data)
        n_nan = eltype(col_data) <: AbstractFloat ? count(isnan, col_data) : 0

        println("Colonne ", col, " : ",
            n_missing, " missing, ",
            n_nothing, " nothing, ",
            n_nan, " NaN")
    end
    #contatenation targetA et B
    #stat.base countmap
    #df .= mapcols(x -> parse.(Float64, x), df)
    plm=PairwiseListMatrices.PairwiseListMatrix(df.rmsd,false,0.0)
    @show first(plm,5)
    mat = PairwiseListMatrices.from_table(df, false, diagonalvalue=1.0)
    @show first(mat,5)
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

function hobohm_clustering(
    table::DataFrames.DataFrame,
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}},
    rmsd_threshold::Float64
)
    rmsds = Dict{Tuple{String,String},Union{Float64,Nothing}}()
    n = DataFrames.nrow(table)
    targets = [table[i, :target] for i in 1:n]
    clustered = Set{String}()
    cluster_labels = Dict{String, Int}()
    cluster_index = 1

    for i in 1:n
        target_i = targets[i]
        if target_i in clustered
            continue
        end

        cluster_labels[target_i] = cluster_index
        push!(clustered, target_i)
        struct_i = structures[target_i]
        if struct_i === nothing
            @warn "Structure for target $target_i is missing."
            continue
        end
        for j in (i+1):n
            target_j = targets[j]
            if target_j in clustered
                continue
            end
            sorted_ids = sort([target_i, target_j])
            key = (sorted_ids[1], sorted_ids[2])

            if !haskey(rmsds, key)
                struct_j = structures[target_j]
                if struct_j === nothing
                    @warn "Structure for target $target_j is missing."
                    continue
                end
                aln = structural_alignment(struct_i, struct_j,
                    BioAlignments.GlobalAlignment(),
                    BioAlignments.AffineGapScoreModel(BioAlignments.BLOSUM62,
                        gap_open=-10, gap_extend=-1))
                if aln !== nothing
                    rmsds[key] = aln[4]
                else
                    rmsds[key] = nothing
                end
            end

            if rmsds[key] !== nothing && rmsds[key] ≤ rmsd_threshold
                cluster_labels[target_j] = cluster_index
                push!(clustered, target_j)
            end
        end

        cluster_index += 1
    end
    @show first(cluster_labels,5)
    # Construction des résultats attendus
    filtered_targets = sort(collect(keys(cluster_labels)))
    clusters = [cluster_labels[t] for t in filtered_targets]
    return cluster_labels, clusters
end

function classical_clustering(
    table::DataFrames.DataFrame,
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}},
    rmsd_threshold::Float64
)

    targets = collect(table.target)
    n = length(targets)

    # Matrice de distances
    dist_matrix = fill(0.0, n, n)

    for i in 1:n
        struct_i = structures[targets[i]]
        if struct_i === nothing
            @warn "Structure for $(targets[i]) missing"
            continue
        end

        for j in i+1:n
            struct_j = structures[targets[j]]
            if struct_j === nothing
                @warn "Structure for $(targets[j]) missing"
                continue
            end

            aln = structural_alignment(
                struct_i, struct_j,
                BioAlignments.GlobalAlignment(),
                BioAlignments.AffineGapScoreModel(
                    BioAlignments.BLOSUM62,
                    gap_open=-10,
                    gap_extend=-1
                )
            )

            if aln !== nothing
                rmsd = aln[4]
                dist_matrix[i,j] = rmsd
                dist_matrix[j,i] = rmsd
            else
                dist_matrix[i,j] = Inf
                dist_matrix[j,i] = Inf
            end
        end
    end

    # Clustering hiérarchique
    hc = Clustering.hclust(dist_matrix, linkage=:average)

    # Découpe selon seuil RMSD
    cluster_ids = Clustering.cutree(hc, h=rmsd_threshold)

    cluster_labels = Dict(targets[i] => cluster_ids[i] for i in 1:n)

    return cluster_labels, cluster_ids
end



function get_target2sequence(expanded_table, msa)
    seqnames=collect(MIToS.MSA.sequencename_iterator(msa))
    clean_seqnames = collect(split(x, '\t')[1] for x in seqnames)
    #=
    seqnames= Set(msa_info)
    @show first(seqnames,5)
    =#
    target2sequence = Dict{String,String}()
    test=[]
    for row in eachrow(expanded_table)
        seqname = ismissing(row.evalue) ? row.query : row.target
        seqname = String(seqname)
        push!(test,seqname)
        if seqname in clean_seqnames
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
        cluster2targets[cluster] = String[]
        for (name, cl) in targets
            if cl == cluster
                if name in cluster2targets[cluster]
                    continue
                end
                push!(cluster2targets[cluster], name)
            end
        end
    end

    cluster2targets, length(cluster2targets)

end

function get_cluster2targets_concatenate(targets, clusters)
    cluster2targets = OrderedCollections.OrderedDict{Int,Vector{String}}()
    # Regroupe les clusters par paquets de 10
    grouped_clusters = Dict{Int, Vector{String}}()
    for (name, cl) in targets
        group = Int(ceil(cl / 10))
        if !haskey(grouped_clusters, group)
            grouped_clusters[group] = String[]
        end
        # Ajoute sans doublon
        if !(name in grouped_clusters[group])
            push!(grouped_clusters[group], name)
        end
    end
    # Convertit en OrderedDict pour compatibilité
    for (group, names) in sort(collect(grouped_clusters))
        cluster2targets[group] = names
    end
    @show length(cluster2targets)
    @show first(cluster2targets, 5)
    cluster2targets, length(cluster2targets)
end

function get_cluster2seqnames(cluster2targets, target2sequence)
    cluster2seqnames = OrderedCollections.OrderedDict{Int,Vector{String}}()
    for (cluster, targets) in cluster2targets
        cluster2seqnames[cluster] = unique(target2sequence[target] for target in targets
                                           if haskey(target2sequence, target))
    end
    @show length(cluster2seqnames)
    cluster2seqnames
end

function _check_names_in_msa(names, msa)
    found_names = []
    seqnames= collect(MIToS.MSA.sequencename_iterator(msa))
    clean_seqnames = collect(split(x, '\t')[1] for x in seqnames)
    for seqname in clean_seqnames
        if seqname in names
            push!(found_names, seqname)
        end
    end
    missing_names = setdiff(names, found_names)
    if !isempty(missing_names)
        @info "The first 5 sequences in the MSA are: $(first(MIToS.MSA.sequencename_iterator(msa), 5))"
        @warn "The following sequences are not in the MSA: $(collect(missing_names))"
    end
    isempty(missing_names)
end

function _rename_msa!(msa)
    new_names = MIToS.MSA.sequencenames(msa)
    MIToS.MSA.rename_sequences!(msa, new_names)
end

function get_cluster2msa(msa, cluster2seqnames)
    query_name = first(MIToS.MSA.sequencename_iterator(msa))
    cluster2msa = OrderedCollections.OrderedDict{Int,MIToS.MSA.AnnotatedMultipleSequenceAlignment}()
    for (cluster, seqnames) in cluster2seqnames
        seqnames = copy(seqnames)
        pushfirst!(seqnames, query_name)
        _rename_msa!(msa)
        if ! _check_names_in_msa(seqnames, msa)
            @warn "The cluster $cluster will be skipped."
            continue
        end
        try
            names=MIToS.MSA.sequencenames(msa)
            all_name= String[]
            for nom in seqnames
                name_line = findfirst(x -> startswith(x, nom), names)
                if name_line !== nothing
                    push!(all_name, names[name_line])
                end
            end
            if isempty(all_name)
                @warn "No MSA for cluster $cluster"
                continue
            end
            @show length(all_name)
            all_name = unique(all_name)
            @show length(all_name)
            cluster2msa[cluster] = msa[all_name, :]
        catch err
            @error "Error ($err) getting the MSA for cluster $cluster"
            @info "seqnames: $seqnames"
            @info "msa names: $(MIToS.MSA.sequencenames(msa))"
            rethrow(err)
        end
    end
    @show length(cluster2msa)
    cluster2msa
end

function get_cluster2structures(structures, cluster2targets)
    # dictionnaire final
    cluster2structures = OrderedCollections.OrderedDict{Int, Dict{String, Vector{MIToS.PDB.PDBResidue}}}()

    # pour tracker les targets déjà assignés
    used_targets = Set{String}()

    for (cluster, targets) in cluster2targets
        # on filtre les targets déjà utilisés
        available_targets = filter(t -> !(t in used_targets), targets)
        
        # on choisit les sous-targets
        subtargets = String[]
        n = length(available_targets)

        if n == 0
            @warn "Cluster $cluster has no available targets left after filtering"
            continue
        elseif n < 8
            subtargets = available_targets[1:min(4, n)]
        else
            step = max(1, div(n, 4))
            for i in 1:4
                idx = i * step
                if idx <= n
                    push!(subtargets, available_targets[idx])
                end
            end
        end

        # marquer ces targets comme utilisés
        foreach(t -> push!(used_targets, t), subtargets)

        # construire cluster2structures
        cluster2structures[cluster] = Dict(
            t => structures[t] for t in subtargets
        )
    end

    @show length(cluster2structures)

    cluster2structures
end


function create_template_clusters_hobohm(
    expanded_table::DataFrames.DataFrame,
    msa::MIToS.MSA.AnnotatedMultipleSequenceAlignment,
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}},
    cutoff::Float64)
    @info "get_target2sequence"
    target2sequence = AlphaConformers.get_target2sequence(expanded_table, msa)
    targets = Set{String}(expanded_table.target)
    @info "cluster_structures"
    #targets, clusters =hobohm_clustering(expanded_table, structures,cutoff)
    output=MIToS.MSA.hobohmI(MIToS.PDB.superimpose, structures, cutoff; threads=true)
    @show output
    return
    @info "get_cluster2targets"
    cluster2targets, nb_cluster = AlphaConformers.get_cluster2targets(targets, clusters)
    #cluster2targets, nb_cluster = AlphaConformers.get_cluster2targets_concatenate(targets, clusters)
    
    @info "get_cluster2seqnames"
    cl2seq = AlphaConformers.get_cluster2seqnames(cluster2targets, target2sequence)
    @info "get_cluster2msa"
    cl2msa = AlphaConformers.get_cluster2msa(msa, cl2seq)
    @info "get_cluster2structures"
    cl2pdb = AlphaConformers.get_cluster2structures(structures, cluster2targets)
    (nb_cluster, cl2msa, cl2pdb)
end


function create_template_clusters(
    expanded_table::DataFrames.DataFrame,
    msa::MIToS.MSA.AnnotatedMultipleSequenceAlignment,
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}},
    cutoff::Float64)
    @info "get_target2sequence"
    target2sequence = AlphaConformers.get_target2sequence(expanded_table, msa)
    targets = Set{String}(expanded_table.target)
    @info "cluster_structures"
    targets, clusters =classical_clustering(expanded_table, structures,cutoff)
    @info "get_cluster2targets"
    cluster2targets, nb_cluster = AlphaConformers.get_cluster2targets(targets, clusters)
    #cluster2targets, nb_cluster = AlphaConformers.get_cluster2targets_concatenate(targets, clusters)
    
    @info "get_cluster2seqnames"
    cl2seq = AlphaConformers.get_cluster2seqnames(cluster2targets, target2sequence)
    @info "get_cluster2msa"
    cl2msa = AlphaConformers.get_cluster2msa(msa, cl2seq)
    @info "get_cluster2structures"
    cl2pdb = AlphaConformers.get_cluster2structures(structures, cluster2targets)
    (nb_cluster, cl2msa, cl2pdb)
end
function prepare_embedding_matrix(df_input, features)

    df = copy(df_input)
    feats = copy(features)

    # log transform evalue
    if :evalue in names(df)
        df.log_evalue = -log10.(df.evalue .+ 1e-300)
    end

    # coverage
    if (:qend in feats) && (:alnlen in feats)
        query_length = maximum(df.qend)
        df.coverage = df.alnlen ./ query_length
        push!(feats, :coverage)
    end

    numeric_features = [f for f in feats if eltype(df[!,f]) <: Number]

    X = Matrix(df[:, numeric_features])

    μ = mean(X, dims=1)
    σ = std(X, dims=1)
    σ[σ .== 0] .= 1

    X = (X .- μ) ./ σ

    return X
end

function umap_emb(df_umap; 
    features=[:query,:target,:fident,:alnlen,:qend,:mismatch,:gapopen,:bits,:evalue,:qtmscore,:ttmscore,:alntmscore,:rmsd], 
    params_umap...)

    X = prepare_embedding_matrix(df_umap, features)

    # normalisation robuste
    μ = mean(X, dims=1)
    σ = std(X, dims=1)

    σ[σ .== 0] .= 1   # éviter division par 0

    X = (X .- μ) ./ σ

    # UMAP
    umap = pyimport("umap")
    embedding = umap.UMAP(;params_umap...).fit_transform(X)

    @info "UMAP embedding completed with $(size(X,1)) samples and $(size(X,2)) features"

    return embedding

end
function plot_dbscan_cluster(embedding)
    
    hdbscan = pyimport("hdbscan")
    @show size(embedding)   # doit être (N, 2)

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=5
    )

    fit = clusterer.fit(embedding)
    labels = fit.labels_
    
    p = plt.scatter(
        embedding[:, 1], embedding[:, 2],embedding[:,3],
        c = :turbo,
        marker_z = labels,     # IMPORTANT
        colorbar = true,
        alpha = 0.25,
        markerstrokewidth = 0,
        markerstrokealpha = 0,
        markersize = 4,
        xlabel = "UMAP 1",
        ylabel = "UMAP 2",
        zlabel = "UMAP 3",
        size = (1000, 800),
        legend = false
    )
    plt.scatter!(p,
        [embedding[1, 1]],
        [embedding[1, 2]],
        [embedding[1, 3]],
        color = :orange, 
        markerstrokewidth = 0,
        markerstrokealpha = 0,
        markersize = 8
    )
    return p, labels
end

function cluster_creation(df_expanded_table,labels)

    unique_clusters = sort(filter(x -> x != -1, unique(labels)))

    println("Clusters found: ", unique_clusters)

    dic_cluster = Dict{Int, Vector{String}}()

    for cid in unique_clusters

        members_idx = findall(labels .== cid)

        println("Cluster $cid size: ", length(members_idx))

        dic_cluster[cid] = String.(df_expanded_table[members_idx, :].target)

    end

    return dic_cluster
end

function create_umap_clusters(
    output_folder::String,
    expanded_table::DataFrame;
    n_neighbors::Int = 15,
    min_dist::Float64 = 0.3
)

    mkpath(output_folder)

    params = (
        n_neighbors = n_neighbors,
        min_dist = min_dist,
        n_components = 3
    )

    embedding = umap_emb(expanded_table; params...)

    plot1, labels = plot_dbscan_cluster(embedding)

    save(joinpath(output_folder, "UMAP_plot.png"), plot1)

    dic_cluster = cluster_creation(expanded_table, labels)

    return dic_cluster
end

function create_folder_structure_hobohm(clusters,
    cl2msa::OrderedCollections.OrderedDict{Int,MIToS.MSA.AnnotatedMultipleSequenceAlignment},
    cl2pdb::OrderedCollections.OrderedDict{Int,Dict{String,Vector{MIToS.PDB.PDBResidue}}};
    out_folder::String=mktempdir())
    #unique_cluster=unique(clusters)
    calpha_template = String[]

    for clust in 1:clusters
        # MSA
        @show clust
        
        cluster_folder = mkdir(joinpath(out_folder, "cluster_$(clust)"))
        msa_file = joinpath(cluster_folder, "sequences.a3m")
        MIToS.MSA.write_file(msa_file, cl2msa[clust], MIToS.MSA.FASTA)

        # Structures
        cluster_template_folder = mkdir(joinpath(cluster_folder, "templates_adaptative"))
        for (target, structure) in cl2pdb[clust]
            #Get the right extension 
            base_name = split(target, ".")[1]
            if startswith(basename(base_name),"AF") 
                name_cif_file=basename(base_name)*".pdb"
                try 
                    MIToS.Utils.download_file("https://alphafold.ebi.ac.uk/files/"*name_cif_file, joinpath(cluster_template_folder,name_cif_file))
                catch e
                    @warn "Problem to download for AFDB the complete template file: "*base_name
                    out_cif = joinpath(cluster_template_folder, name_cif_file)
                    MIToS.PDB.write_file(out_cif, structure, MIToS.PDB.PDBFile)
                    push!(calpha_template, out_cif)

                    continue
                end
                
            else 
                pdbcode = split(basename(base_name), '.')[1]
                pdbcode = split(basename(pdbcode), '_')[1]
                chain_id = lowercase(split(basename(base_name), '_')[end])
                @show pdbcode
                name_cif_file=pdbcode*".cif"
                try 
                    MIToS.Utils.download_file("https://files.rcsb.org/download/"*name_cif_file, joinpath(cluster_template_folder,name_cif_file))

                catch e
                    @warn "Problem to download from PDB the template file: "*base_name
                    out_cif = joinpath(cluster_template_folder, name_cif_file)
                    MIToS.PDB.write_file(out_cif, structure, MIToS.PDB.MMCIFFile)
                    push!(calpha_template, out_cif)
                    continue
                end
                
            end
        end
        
    end
    CSV.write(joinpath(out_folder, "calpha_template.csv"),
              (; file = calpha_template))
    out_folder
end

function write_new_msa(new_msa_path,query_id,query_seq,cluster_ids,cluster_seqs)
    open(new_msa_path, "w") do io
        println(io, ">$query_id")
        println(io, query_seq)
        for (id, seq) in zip(cluster_ids, cluster_seqs)
            println(io, ">$id")
            println(io, seq)
        end
    end
end


function get_new_template(cluster_ids)
    # on choisit les sous-targets
    templates = String[]
    n = length(cluster_ids)

    if n == 0
        @warn "Cluster $cluster has no available targets left after filtering"
        return nothing
    elseif n < 8
        templates = cluster_ids[1:min(4, n)]
    else
        step = max(1, div(n, 4))
        for i in 1:4
            idx = i * step
            if idx <= n
                push!(templates, cluster_ids[idx])
            end
        end
    end

    return templates
end


function download_template(cluster_folder,templates)
    cluster_template_folder = mkdir(joinpath(cluster_folder, "templates_complete"))
    for target in templates
        #Get the right extension 
        base_name = split(target, ".")[1]
        if startswith(basename(base_name),"AF") 
            new_name= replace(basename(base_name), "MODEL_V6" => "model_v6")
            name_cif_file=basename(new_name)*".pdb"
            try 
                MIToS.Utils.download_file("https://alphafold.ebi.ac.uk/files/"*name_cif_file, joinpath(cluster_template_folder,name_cif_file))
            catch e
                @warn "Problem to download for AFDB the complete template file: "*name_cif_file
                continue
            end
        else 
            pdbcode = split(basename(base_name), '.')[1]
            pdbcode = split(basename(pdbcode), '_')[1]
            chain_id = lowercase(split(basename(base_name), '_')[end])
            @show pdbcode
            name_cif_file=pdbcode*".cif"
            try 
                MIToS.Utils.download_file("https://files.rcsb.org/download/"*name_cif_file, joinpath(cluster_template_folder,name_cif_file))

            catch e
                @warn "Problem to download from PDB the template file: "*name_cif_file
                continue
            end
            
        end
    end
    
end

function get_files_cluster(new_cluster,query_id,query_seq,cluster_ids,cluster_seqs)
    new_msa_path=joinpath(new_cluster,"sequences.a3m")
    write_new_msa(new_msa_path,query_id,query_seq,cluster_ids,cluster_seqs)

    #Create template folder
    template_names=get_new_template(cluster_ids)
    download_template(new_cluster,template_names)
end


function chunk_sequences(seqs,
                         ids;
                         min_size::Int = 5,
                         max_size::Int = 10)

    
    n = length(seqs)
    if n == 0
        return Vector{Vector{String}}(), Vector{Vector{String}}()
    end


    # nombre minimal de clusters pour ne pas dépasser max_size
    n_clusters = ceil(Int, n / max_size)

    # taille de base de chaque cluster
    base_size = div(n, n_clusters)

    # combien de clusters auront +1 élément
    remainder = n % n_clusters

    cluster_seqs = Vector{Vector{String}}()
    cluster_ids  = Vector{Vector{String}}()

    i = 1

    for c in 1:n_clusters

        size = base_size + (c <= remainder ? 1 : 0)

        push!(cluster_seqs, seqs[i:i+size-1])
        push!(cluster_ids,  ids[i:i+size-1])

        i += size
    end

    return cluster_seqs, cluster_ids
end

function create_folder_structure_subdivide(clusters,main_new_folder)
    all_seq_a3m_path=joinpath(main_new_folder, "all_sequence.a3m")
    ids,sequences =read_a3m(all_seq_a3m_path)
    query_id=ids[1]
    query_seq=sequences[1]
    cluster_id=1
    check_files=[]
    for (cid, members) in clusters
        @show cid
        @show members
        if length(members) >= 4 && length(members) <= 9
            
            #If have the right number of sequences => keep cluster as it is 
            new_cluster=joinpath(main_new_folder,"cluster_$cluster_id")
            mkdir(new_cluster)
            idx = findall(id -> id in members, ids)

            cluster_seqs = sequences[idx]
            @show length(cluster_seqs)
            @show length(members)
            get_files_cluster(new_cluster,query_id,query_seq,members,cluster_seqs)
            cluster_id=cluster_id+1

        elseif length(members)>9
            #If have more sequence => divide to be between 5 and 10
            idx = findall(id -> id in members, ids)

            cluster_all_seqs = sequences[idx]
            cluster_seqs,cluster_ids=chunk_sequences(cluster_all_seqs,members,min_size=4,max_size=9)
            
            #recreate the cluster
            for i in 1:length(cluster_ids)
                new_cluster=joinpath(main_new_folder,"cluster_$cluster_id")
                mkdir(new_cluster)
                get_files_cluster(new_cluster,query_id,query_seq,cluster_ids[i],cluster_seqs[i])
                cluster_id=cluster_id+1
            end
        else 
            append!(check_files,members)
        end
    end
    idx_outlier = findall(id -> id in check_files, ids)
    cluster_all_seqs = sequences[idx_outlier]
    cluster_seqs,cluster_ids=chunk_sequences(cluster_all_seqs,check_files,min_size=4,max_size=9)
    for i in 1:length(cluster_ids)
        new_cluster=joinpath(main_new_folder,"cluster_$cluster_id")
        mkdir(new_cluster)
        get_files_cluster(new_cluster,query_id,query_seq,cluster_ids[i],cluster_seqs[i])
        cluster_id=cluster_id+1
    end

    return cluster_id-1
end 


function create_folder_structure(large_small_pairs::Vector{Tuple{Int,Int}},
    large_cl2msa::OrderedCollections.OrderedDict{Int,MIToS.MSA.AnnotatedMultipleSequenceAlignment},
    small_cl2pdb::OrderedCollections.OrderedDict{Int,Dict{String,Vector{MIToS.PDB.PDBResidue}}};
    out_folder::String=mktempdir())
    for (large, small) in large_small_pairs
        # MSA
        cluster_folder = mkdir(joinpath(out_folder, "cluster_$(large)_$(small)"))
        msa_file = joinpath(cluster_folder, "sequences.a3m")
        MIToS.MSA.write_file(msa_file, large_cl2msa[large], MIToS.MSA.FASTA)
        # Structures
        cluster_template_folder = mkdir(joinpath(cluster_folder, "templates"))
        for (target, structure) in small_cl2pdb[small]
            #Get the right extension 
            base_name = replace(target, ".pdb" => "")
            if startswith(base_name, "AF")
                upper_file = joinpath(cluster_template_folder, uppercase(base_name) * ".pdb")
                MIToS.PDB.write_file(upper_file, structure, MIToS.PDB.PDBFile)
            else 
                # lowercase and uppercase file 
                lower_file = joinpath(cluster_template_folder, lowercase(base_name) * ".pdb")
                MIToS.PDB.write_file(lower_file, structure, MIToS.PDB.PDBFile)
            end
        end
    end
    out_folder
end

#Clean the file name in the MSA and prep the template file for AlphaFold (.cif + good name)
function clean_msa_template_names(clusters,out_folder;cluster_name::String="templates_adaptative")

    for clust in 1:clusters
        corresponding_name=DataFrames.DataFrame(Ref_name=String[],New_name_msa=String[],New_name_template=String[]) # Create an empty DF

        # MSA
        @show clust

        cluster_folder = joinpath(out_folder, "cluster_$(clust)")
        msa_path = joinpath(cluster_folder, "sequences.a3m")

        #1- read the msa
        ids, sequences = read_a3m(msa_path)
        query_id = ids[1]
        query_seq = sequences[1]
        
        #2 - Clean name
        clean_ids=String[]
        for name in ids
            if name == query_id
                push!(clean_ids, query_id)
            elseif startswith(basename(name), "AF")
                push!(clean_ids, uppercase(name))
            else
                pdb_id = split(basename(name), ".")[1]
                chain_id = split(basename(name), "_")
                if length(chain_id) > 1
                    chain_id = chain_id[end]
                    
                else
                    chain_id = "a"
                end
                push!(clean_ids, lowercase("$(pdb_id)_$(chain_id)"))
            end
        end
        
        #=
        #3- Rewrite MSA with clean names
        open(msa_path, "w") do io
            for (id, seq) in zip(clean_ids, sequences)
                println(io, ">$id")
                println(io, seq)
            end
        end
        =#
        cluster_template_folder = joinpath(cluster_folder, "templates_complete")
        if isdir(joinpath(cluster_folder, cluster_name))
            @info "Folder already clean"
            continue
        end
        cluster_template_clean_folder = mkdir(joinpath(cluster_folder, cluster_name))
        
        templates = glob("*",cluster_template_folder)
        isempty(templates) && error("No templates found in $cluster_template_folder")
        i=0
        for template in templates 
            i+=1
            #4- Change template name 
            if startswith(basename(template), "AF")
                clean_template_name = "t00$(i)_a"
                file_name="t00$(i).pdb"
                push!(corresponding_name,(basename(template),clean_template_name,file_name))
                name=split(basename(template), ".")[1]
                for i in 1:length(clean_ids)
                    if clean_ids[i]== uppercase(name)
                        clean_ids[i] = clean_template_name
                    end
                end
            else
                clean_template_name = basename(template)
                checks = split(clean_template_name, "_")
                name=split(clean_template_name, "_")[1]
                if length(checks) < 2
                    clean_template_name = name*"_a"
                end
                base_name = replace(name, ".cif" => "")
                file_name=lowercase(base_name)*".cif"
            end
            cp(template, joinpath(cluster_template_clean_folder,file_name); force=true)
        end    

        
        #5- Rewrite MSA with clean names
        open(msa_path, "w") do io
            for (id, seq) in zip(clean_ids, sequences)
                println(io, ">$id")
                println(io, seq)
            end
        end
            
        CSV.write(cluster_folder*"/corresponding_name_template.csv", corresponding_name) 
    end
end

function write_a3m(filename::String, mafft_msa::Vector{Any})
    open(filename, "w") do io
        for ( header,seq) in mafft_msa
            println(io, ">$header")
            println(io, seq)
        end
    end
end

function align_mafft(msa_file::String, pdb_dir::String)
    #Get the MSA file output by Foldseek
    msas = MIToS.MSA.read_file(msa_file, MIToS.MSA.FASTA, generatemapping=true)
    seqnames = MIToS.MSA.sequencenames(msas) # Get the sequence names from the MSA

    mafft_msa = []
    model = "1" # Default model
    query= first(seqnames) # The first sequence is the query sequence
    pdb_id =String(split(query, "_")[1])
    chain_id =String(split(query, "_")[2])

    struc = MIToS.PDB.read_file(joinpath(pdb_dir,pdb_id)*".pdb", MIToS.PDB.PDBFile)
    @show typeof(struc)
    sequences=MIToS.PDB.modelled_sequences(struc;chain=chain_id)
    seq = sequences[(model=model, chain=chain_id)]

    push!(mafft_msa, (query, seq))
    for seqname in seqnames[2:end]

        pdb_id =split(seqname, ".")[1] # Get the PDB file name from the sequence name
        struc = MIToS.PDB.read_file(joinpath(pdb_dir,pdb_id)*".pdb", MIToS.PDB.PDBFile)

        chain_id = split(seqname, "_") # Get the chain from the sequence name
        if length(chain_id) > 1
            chain_id = String(split(chain_id[2],"\t")[1]) # Get the chain ID from the sequence name
            sequences=MIToS.PDB.modelled_sequences(struc;chain=chain_id)
            seq = sequences[(model=model, chain=chain_id)]
        else
            sequences=MIToS.PDB.modelled_sequences(struc)
            seq = first(values(sequences))
        end
        
        push!(mafft_msa, (seqname, seq))
    end
    @show length(mafft_msa)
    @show first(mafft_msa, 5)
    @show typeof(mafft_msa)
    write_a3m("tmp_for_mafft.fasta", mafft_msa)
    # Align the MSA using MAFFT

    @info "Aligning MSA with MAFFT: $msa_file"
    
    run(pipeline(`$(MAFFT_jll.mafft()) --auto tmp_for_mafft.fasta`, stdout=msa_file))
    @show typeof(msa_file)
    @show msa_file
    # Lire l'alignement aligné
    aligned_msa = MIToS.MSA.read_file(msa_file, MIToS.MSA.FASTA)
    # Adjust the reference sequence in the MSA to avoid gaps
    MIToS.MSA.adjustreference!(aligned_msa)
    # Write the aligned MSA to the output file
    MIToS.MSA.write_file(msa_file, aligned_msa, MIToS.MSA.A3M)
end

function get_res_beg_end(pdb_file::String)
    res_ids = Int[]

    open(pdb_file, "r") do io
        for line in eachline(io)
            if startswith(line, "ATOM")
                raw = strip(line[23:26])

                # 👉 garder uniquement les chiffres
                num = match(r"\d+", raw)

                if num !== nothing
                    push!(res_ids, parse(Int, num.match))
                end
            end
        end
    end

    return minimum(res_ids), maximum(res_ids)
end

function align_full_seq(full_seq::String, msa::Vector{String})
    ref_seq = msa[1]

    aligned_ref = ""
    mapping = Int[]  # positions dans ref_seq ou 0 si gap

    i = 1  # index ref_seq

    for c in full_seq
        if i <= length(ref_seq) && c == ref_seq[i]
            push!(mapping, i)
            i += 1
        else
            push!(mapping, 0)  # gap
        end
    end

    # Construire nouveau MSA
    new_msa = String[]

    # 1️⃣ première ligne = full_seq (ta référence propre)
    push!(new_msa, full_seq)

    # 2️⃣ autres séquences
    for seq in msa[2:end]
        new_seq = ""
        for pos in mapping
            if pos == 0
                new_seq *= "-"
            else
                new_seq *= seq[pos]
            end
        end
        push!(new_msa, new_seq)
    end

    return new_msa
end

#faire pour ajouter n pdb et les concaténer dans foldseek 
function alphaconformers(input_pdb, pdb_folder, out_folder,n_threads; db::Vector{String}=["/alpha/database/pdb/fullpdb"], 
        evalue_cutoff::Float64=1e-5,cutoff::Float64=1.0, mafft::Bool=false,
        umap::Bool=false,n_neighbors::Int = 15,min_dist::Float64 = 0.3,
        test_analyse::Bool=false,
        full_seq::String=nothing
        )
    @info "Running Foldseek"
    query_pdb_code = split(basename(input_pdb),"_")[1]
    query_chain_code = split(split(basename(input_pdb),"_")[2],".")[1]
    
    output = run_foldseek(input_pdb, n_threads, db, out_folder=out_folder) 
    # Initialiser un vecteur vide de String
    output_vector = Vector{String}()
    
    # Boucle correcte sur les éléments de `output`
    for item in output
        if mafft ## Align MSA on sequence and not on structure 
            align_mafft(item.msa_file,pdb_folder)  # Aligner les séquences MSA
        end 
        push!(output_vector, item.table_file)  # Ajouter au vecteur
    end
    
    # Fusionner les tables
    merged_table = merge_tables(output_vector)
    @show size(merged_table)
    if !isnan(evalue_cutoff)
        filter!(row -> row.evalue < evalue_cutoff, merged_table)
    end
    
    @show size(merged_table)
    if size(merged_table, 1) == 0
        @error "No results found after filtering with e-value cutoff: $evalue_cutoff"
        return
    end
    @info "Adding known conformations"
    @show size(merged_table)
    sifts_uniprot_mapping = get_uniprot_mapping()
    clean_table = add_known_conformations!(deepcopy(merged_table), sifts_uniprot_mapping,pdb_folder,out_folder,input_pdb,n_threads)
    @show size(clean_table)
    
    if test_analyse
        expanded_table=delete_query_from_target(deepcopy(clean_table),sifts_uniprot_mapping,String(query_pdb_code),String(query_chain_code))
        @show size(expanded_table)
        isempty(expanded_table) && error("No targets remaining after removing query from targets.")
        CSV.write(joinpath(out_folder,"expanded_table_delete.csv"),expanded_table)
    else 
        expanded_table=deepcopy(clean_table)
    end
    
    @info "Merge MSAS"
    merged_msa = merge_msas(expanded_table)
    MIToS.MSA.write_file(joinpath(out_folder, "all_sequence.a3m"), merged_msa, MIToS.MSA.A3M)
    @show size(merged_msa)

    ### Align sequence on AFDB sequences - to not have missing residues 
    if full_seq !== nothing
        res_beg,res_end=get_res_beg_end(input_pdb)
        @show res_beg,res_end
        same_sull_seq=full_seq[res_beg:res_end]
        @show length(full_seq)
        @show length(same_sull_seq)
        align_merged_msa=align_full_seq(same_sull_seq,merged_msa)
        MIToS.MSA.write_file(joinpath(out_folder, "all_sequence_align_uniprot_seq.a3m"), align_merged_msa, MIToS.MSA.A3M)
    else
        align_merged_msa=deepcopy(merged_msa)
    end
    return
    if !umap
        #merged_msa = merge_msas(merged_table)
        @info "Getting the aligned structures"
        structures = _get_aligned_structures(expanded_table)
        @show length(structures)

        
        #uniprot2targets = get_uniprot2targets(target2uniprot, expanded_table)
        @show size(expanded_table)
        
        @info "Measuring RMSDs"
        #rmsds = process_known_conformations!(structures, expanded_table, target2uniprot, uniprot2targets, pdb_folder=pdb_folder)
        #@show size(expanded_table)
        #fill_rmsds!(rmsds, expanded_table, structures)
        
        
        @info "Clustering structures"
        clusters, cl2msa, cl2pdb= create_template_clusters_hobohm(expanded_table, align_merged_msa, structures,cutoff)
        #clusters, cl2msa, cl2pdb = create_template_clusters(expanded_table, merged_msa, structures,cutoff)
        @info "Create folder"
        #create_folder_structure(large_small_pairs, large_cl2msa, small_cl2pdb, out_folder=out_folder)
        
        create_folder_structure_hobohm(clusters, cl2msa, cl2pdb, out_folder=out_folder)
    else 

        @info "Clustering with UMAP"
        dic_clusters=create_umap_clusters(out_folder,expanded_table;n_neighbors,min_dist)
        clusters=create_folder_structure_subdivide(dic_clusters,out_folder)
    end

    clean_msa_template_names(clusters,out_folder)
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