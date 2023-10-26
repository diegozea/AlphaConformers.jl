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

"""
    read_foldseek_search_results(file::AbstractString)

Reads the Foldseek easy-search output file (m8) and returns a DataFrame with the results 
and proper column names.
"""
function read_foldseek_search_results(file::AbstractString)
    DataFrames.DataFrame(CSV.File(file, delim='\t', header=_M8_COL_NAMES))
end

function _create_folder(folder)
    if isdir(folder)
        @warn "$folder already exists. Results will be overwritten."
        rm(folder; recursive=true)
    end
    mkdir(folder)
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
    _create_folder(out_folder_db)
    pdb_name = first(splitext(basename(pdb_file))) # filename without extension
    table_file = joinpath(out_folder_db, "$(pdb_name)_results.m8")
    msa_file = joinpath(out_folder_db, "msa.a3m")
    aligned_structures_folder = joinpath(out_folder_db, "aligned_structures")
    _create_folder(aligned_structures_folder)
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
