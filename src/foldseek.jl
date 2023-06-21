# Functions to run Foldseek inside the pipeline

ENV["FOLDSEEK_DB_PATH"] = ""

"""
    foldseek_search(pdb_file::AbstractString; db_path::String = ENV["FOLDSEEK_DB_PATH"])

Searches a given protein structure in the Foldseek database. 

The `pdb_file` is the path to the protein structure file in PDB format. The `db_path` is 
the path to the Foldseek database, which defaults to the `FOLDSEEK_DB_PATH` environment 
variable. If the `FOLDSEEK_DB_PATH` environment variable is not set and `db_path` is not 
given, the function will throw an error. 

The function executes a Foldseek easy-search, saving the results in an output file with 
the same base name as the PDB file but with `_results.m8` appended to it. The function 
operates in a temporary directory that is automatically cleaned up afterwards. 

The function returns the path to the Folseek easy-search output file.
"""
function foldseek_search(pdb_file::AbstractString; db_path::String = ENV["FOLDSEEK_DB_PATH"])
    isempty(db_path) && error("Please set the FOLDSEEK_DB_PATH environment variable or " * 
        "the db_path keyword argument to the path of the Foldseek database.")
    isfile(db_path) || error("The path to the Foldseek database is not a file.")
    mktempdir() do tmp_folder
        # path without extension
        path = first(splitext(abspath(pdb_file)))
        out_file = "$(path)_results.m8"
        run(`$(Foldseek_jll.foldseek()) easy-search $pdb_file $db_path $out_file $tmp_folder`)
        out_file
    end
end

