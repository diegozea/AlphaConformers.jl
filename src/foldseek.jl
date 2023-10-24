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

