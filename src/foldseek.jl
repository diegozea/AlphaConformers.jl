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
        out_file = "$(path)_results.m8"
        run(`$(Foldseek_jll.foldseek()) easy-search $pdb_file $db_path $out_file $tmp_folder --format-output $format_output`)
        
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
    else
        out_folders = table
    end
    
    msas = [MIToS.MSA.read_file(joinpath(folder, "msa.a3m"), MIToS.MSA.FASTA, generatemapping=true)
            for folder in out_folders]
    # Select only the matched sequences by using bits and fident to identify the matches.
    # Apply this filter only when there are duplicated names to prevent losing elements 
    # due to numerical differences in comparisons.
    if size(msas)[1] == 0
        @error "No MSAs found in the provided table."
        return nothing
    end
    for i in eachindex(msas)
        msa = msas[i]
        
        seqnames = MIToS.MSA.sequencenames(msa)[2:end]
        
        if abs(length(seqnames) - size(msa, 1)) != 1
            @warn "Warning, the number of sequence names does not match the number of rows in the MSA!"
            while length(seqnames) < (size(msa, 1)-1)
                push!(seqnames, "Unnamed_" * string(length(seqnames) + 1))
            end
        end
        
        msa_targets = [first(split(seqname, "\t")) for seqname in seqnames]
        duplicated_msa_targets = _find_duplicates(msa_targets)
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
        msas[i] = msa[selected, :]
    end
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
    # Clean the sequence names by deleting the MSA number at the beginning that 
    # was added by the join function and the sequence data Foldseek adds at the end.
    msa_a_seq=collect(MIToS.MSA.sequencename_iterator(msa_a))
    cleaned_names=String[]
    for name in msa_a_seq
        all_name=split(name,"\t")[1]
        push!(cleaned_names, join(split(all_name, "_")[2:end], "_"))
    end 
    
    if length(cleaned_names) != length(unique(cleaned_names))
        @info "Duplicate names found in the MSA. Removing duplicates by keeping the one with the best score."

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









