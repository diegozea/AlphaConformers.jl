

function _read_pdb(file::String, chain_code)
    # Read the chain specified by chain_code from the PDB file
    # occupancyfilter is needed to avoid the duplicated residue warnings with TMalign
    read(file, MIToS.PDB.PDBFile, chain=chain_code,
        onlyheavy=true, occupancyfilter=true)
end


function _get_pdb_chain(pdb_db::String, pdb_code, chain_code; extention=".pdb")
    # If the pdb_db path is an empty string, download the PDB file from the web
    chain = nothing
    if isempty(pdb_db)
        # Download the PDB file from the web into a temporary file
        tmp_path = joinpath(tempdir(), pdb_code * ".pdb.gz")
        try
            MIToS.PDB.downloadpdb(pdb_code, filename=tmp_path, format=MIToS.PDB.PDBFile)
            chain = _read_pdb(tmp_path, chain_code)
        catch err
            @error "Error downloading or reading the PDB file $pdb_code: $err"
        finally
            isfile(tmp_path) && rm(tmp_path) # Delete the temporary file
        end
    else
        pdb_path = joinpath(pdb_db, pdb_code * extention)
        if isfile(pdb_path)
            chain = _read_pdb(pdb_path, chain_code)
        end
    end
    chain
end

"""
    _create_empty_folder(path)

Helper function to create an empty folder at the given path. If the folder already exists,
all its contents will be deleted.
"""
function _create_empty_folder(path)
    if isdir(path)
        @warn "The folder $path already exists; all its contents will be deleted."
        rm(path, recursive=true, force=true)
    end
    @info "Creating folder $path"
    mkdir(path)
end

"""
    create_pdb_folder(targets::Union{Set, DataFrames.DataFrame, AbstractArray}, 
                      folder_path::String; 
                      pdb_db::String = get(ENV, "PDB_DB", ""), 
                      extension::String = ".pdb")

Creates a local folder at `folder_path` and populates it with the PDB files specified by 
the `targets`. The `targets` argument can be a Set, an AbstractArray, or a DataFrame with 
a `target` column, such as the results loaded from a Foldseek search using the 
`read_foldseek_search_results` function. Each target represents a PDB code, possibly 
followed by a chain code. If a chain code is specified, only that chain is saved in the 
new folder. If no chain code is specified, all chains are saved. The expected format of 
the `targets` depends on the number of chains:
  * 1 chain:    "pdb"
  * 2 chains:   "pdb_chain"
For example: "1EX7.pdb" or "4F4J.pdb_B". Note that the presence of the "pdb" extension is 
not mandatory and depends on the database used when reading Foldseek results.

The PDB files are either fetched from a local database path specified by `pdb_db` or 
downloaded from the web if `pdb_db` is an empty string. Note that you can set the `pdb_db` 
keyword argument when calling the function, or alternatively, set the `PDB_DB` 
environment variable. When retrieving the PDB files from a local folder, the optional 
`extension` argument specifies the file extension of the PDB files in that folder 
(default is ".pdb").

!!! warning
    If the specified `folder_path` already exists, it will be deleted and a new folder 
    will be created.

The function returns the absolute path to the created folder.
"""
function create_pdb_folder(targets::Set, folder_path::String; 
        pdb_db::String = get(ENV, "PDB_DB", ""), extention::String = ".pdb")
    local_pdb_folder = abspath(folder_path)

    # If the folder already exists, delete it and create a new one
    _create_empty_folder(folder_path)
    
    for target in targets
        # Split target into pdb_code and chain_code.
        # NOTE: Foldseek outputs targets in different formats based on the chain count:
        #   1 chain:    pdb (e.g. 1EX7.pdb)
        #   2 chains:   pdb_chain (e.g. 4F4J.pdb_B)
        # NOTE: The presence of the 'pdb' extension depends on the Foldseek database used.
        fields = split(String(target), '_')
        pdb_code = first(splitext(String(fields[1]))) # remove the extension
        @assert length(pdb_code) == 4 "The PDB code must be 4 characters long."
        if length(fields) == 2
            chain_code = String(fields[2])
        else
            chain_code = MIToS.PDB.All
        end
        # Define the file path to the PDB within the newly created folder.
        # The file name is the same as the one in the list of targets.
        local_pdb_path = joinpath(local_pdb_folder, String(target))
        # Read the original PDB file and save the chain in the new folder.
        chain = _get_pdb_chain(pdb_db, pdb_code, chain_code; extention=extention)
        if chain !== nothing
            MIToS.PDB.write(local_pdb_path, chain, MIToS.PDB.PDBFile)
        end
    end

    return local_pdb_folder
end

function create_pdb_folder(targets::DataFrames.DataFrame, folder_path::String; 
        pdb_db::String = get(ENV, "PDB_DB", ""), extention::String = ".pdb")
    create_pdb_folder(Set{String}(targets.target), folder_path; 
        pdb_db=pdb_db, extention=extention)
end

function create_pdb_folder(targets::AbstractArray, folder_path::String; 
        pdb_db::String = get(ENV, "PDB_DB", ""), extention::String = ".pdb")
    create_pdb_folder(Set{String}(targets), folder_path; pdb_db=pdb_db, extention=extention)
end