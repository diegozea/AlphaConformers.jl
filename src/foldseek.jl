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
        db_path::String = get(ENV, "FOLDSEEK_DB_PATH", ""),
        format_mode::Int = 0)
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
        db_path::String = get(ENV, "FOLDSEEK_DB_PATH", "");
        out_folder::String = dirname(abspath(pdb_file)))
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

    [(;table_file, msa_file, aligned_structures_folder)]
end

function run_foldseek(pdb_file::AbstractString, db_path::Vector{String};
    out_folder::String = dirname(abspath(pdb_file)))
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
    (
        table = table_file,
        msa = joinpath(table_folder, "msa.a3m"),
        structures = joinpath(table_folder, "aligned_structures")
    )
end

# concatenate msas ----------------------------------------------------------------------- #

function merge_msas(table::DataFrames.DataFrame)
    out_folders = dirname.(unique(table.table_file))
    msas = [ read(joinpath(folder, "msa.a3m"), MIToS.MSA.FASTA) 
        for folder in out_folders ]
    # Return the MSA if there is only one
    length(msas) == 1 && return only(msas)
    # Otherwise, concatenate the multiple MSAs using the query sequence as reference
    nothing
end

function _match_columns(msa_a, msa_b, ref_a, ref_b)
    ref_a = replace(MIToS.MSA.stringsequence(msa_a, ref_a), "-" => "")
    ref_b = replace(MIToS.MSA.stringsequence(msa_b, ref_b), "-" => "")
    aln_model = BioAlignments.AffineGapScoreModel(match=6, mismatch=-4, 
        gap_open=-2, gap_extend=-1)
    aln = BioAlignments.pairalign(BioAlignments.GlobalAlignment(), ref_a, ref_b, aln_model)
    pos_a = 0
    pos_b = 0
    positions_a = Int[]
    positions_b = Int[]
    for (res_a, res_b) in alignment(aln)
        if res_a != '-'
            pos_a += 1
            push!(positions_a, pos_a)
        else
            push!(positions_a, 0)
        end
        if res_b != '-'
            pos_b += 1
            push!(positions_b, pos_b)
        else
            push!(positions_b, 0)
        end
    end
    (positions_a, positions_b)
end

function _define_seq_blocks(positions_a, positions_b)
    blocks = Tuple{Vector{Int}, Vector{Int}}[(Int[], Int[]),]
    previous_a = positions_a[1]
    previous_b = positions_b[1]
    for (a, b) in zip(positions_a, positions_b)
        if ((a == 0 && previous_a != 0) || 
            (b == 0 && previous_b != 0) || 
            (a != 0 && b != 0 && (previous_a == 0 || previous_b == 0)))
            push!(blocks, (Int[a], Int[b]))
        else
            push!(blocks[end][1], a)
            push!(blocks[end][2], b)
        end
        previous_a, previous_b = a, b
    end
    blocks
end





# structure_distances -------------------------------------------------------------------- #


function _aligned_pdb_file(query, target, aligned_structures_folder)
    file_name = "aln_$(query)_$(target).pdb"
    joinpath(aligned_structures_folder, file_name)
end



function get_aligned_residues(query, target, paths)
    msa = read(paths.msa, MIToS.MSA.FASTA)
    
    pdb_file = _aligned_pdb_file(query, target, paths.structures)
    res = read(pdb_file, MIToS.PDB.PDBFile)

end


function structure_distances(foldseek_results::DataFrames.DataFrame)
    targets = foldseek_results.target
    n_targets = length(targets) # also nrows
    distances = zeros(Float16, n_targets, n_targets)
    for i in 1:(n_targets-1)
        row_i = foldseek_results[i, :]
        paths_i = _get_paths(row_i.table_file)
        pdb_file_i = _aligned_pdb_file(row_i.query, row_i.target, paths_i.structures)
        for j in (i+1):n_targets
            row_j = foldseek_results[j, :]
            paths_j = _get_paths(row_j.table_file)
            pdb_file_j = _aligned_pdb_file(row_j.query, row_j.target, paths_j.structures)
            distances[i, j] = distances[j, i] = MIToS.PDB.rmsd(pdb_a, pdb_b)
        end
    end
end



