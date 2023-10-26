# This file defines a series of utils to use the UniProt annotations from SIFT to 
# identify known conformations from a series of PDB files. Using the UniProt mapping 
# also allows to filter out structures from the query protein to evaluate the pipeline.

# UniProt mapping
# ===============

"""
    get_uniprot_mapping(sifts_db::String = get(ENV, "SIFTS_DB", pwd()))

Reads the `pdb_chain_uniprot.csv.gz` file and returns a DataFrame. If the file is not found 
in the specified `sifts_db` folder, it is downloaded from the SIFTS FTP server. 
When calling the function, you can set the `sifts_db` positional argument or the 
`SIFTS_DB` environment variable. If this is not done, the current working directory will be
used by default.
"""
function get_uniprot_mapping(sifts_db::String = get(ENV, "SIFTS_DB", pwd()))
    @assert isdir(sifts_db) "sifts_db must be a folder containing pdb_chain_uniprot.csv.gz"
    sifts_file_name = "pdb_chain_uniprot.csv.gz"
    sifts_file_path = joinpath(sifts_db, sifts_file_name)
    if !isfile(sifts_file_path)
        @info "Downloading pdb_chain_uniprot.csv.gz into $sifts_db"
        url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz"
        Downloads.download(url, sifts_file_path)
    end
    DataFrames.DataFrame(CSV.File(sifts_file_path, 
        comment="#", missingstring=["", "None"]))
end

"""
    get_uniprot_acc(data::DataFrames.DataFrame, 
        pdb::String, chain::Union{String, Type{MIToS.PDB.All}}=MIToS.PDB.All)

This function returns a vector of UniProt accessions associated with a given PDB code and
chain (if available). If the `chain` argument is not specified or if it is set to
`MIToS.PDB.All` all UniProt accessions associated with the PDB code will be returned.

This function takes as input the `DataFrame` returned by the `get_uniprot_mapping` 
function as its first positional argument, `data`.
"""
function get_uniprot_acc(data::DataFrames.DataFrame, 
        pdb::String, chain::Union{String, Type{MIToS.PDB.All}}=MIToS.PDB.All)
    pdb_code = lowercase(pdb)
    if chain !== MIToS.PDB.All
        ups = data[(data.PDB .== pdb_code) .& (data.CHAIN .== chain), :SP_PRIMARY]
    else
        ups = data[data.PDB .== pdb_code, :SP_PRIMARY]
    end
    String.(unique(ups))
end

"""
    get_pdb_codes(data::DataFrames.DataFrame, uni_acc::String)

This function returns a `DataFrame` with two columns, `PDB` and `CHAIN`, containing the
PDB codes and chains associated with a given UniProt accession as uppercase strings.
This function takes as input the `DataFrame` returned by the `get_uniprot_mapping` 
function as its first positional argument, `data`.
"""
function get_pdb_codes(data::DataFrames.DataFrame, uni_acc::String)
    pdb_codes = data[data.SP_PRIMARY .== uni_acc, [:PDB, :CHAIN]]
    uppercase.(unique(pdb_codes))
end

# FoldSeek results
# ================

# Functions to filter out the query protein from the FoldSeek search results

function _get_pdb_and_chain(pdb_file::AbstractString)
    pdb_chain = replace(uppercase(basename(pdb_file)), 
        ".PDB" => "", ".GZ" => "", 
        # delete anythig that is not an uppercase letter or a number
        r"[^A-Z0-9]" => "")
    # We expect the PDB code to be the first 4 characters
    if length(pdb_chain) > 4
        chain = String(pdb_chain[5:end])
    else
        chain = MIToS.PDB.All
    end
    String(pdb_chain[1:4]), chain
end

function _is_chain(pdb_file::String, pdb_code::String, chain_code::String)
    # Actual pdb and chain codes determined by the pdb_file
    actual_pdb, actual_chain = _get_pdb_and_chain(pdb_file)
    if actual_chain === MIToS.PDB.All
        actual = uppercase(string(actual_pdb))
        expected = uppercase(string(pdb_code))
    else
        actual = uppercase(string(actual_pdb, actual_chain))
        expected = uppercase(string(pdb_code, chain_code))
    end
    actual == expected
end

# delete the query if it is in the target column
function delete_query_from_target!(search_results::DataFrames.DataFrame, 
        sifts_uniprot_mapping::DataFrames.DataFrame,
        query_pdb_code::String, query_chain_code::String)
    query_uniprots = get_uniprot_acc(sifts_uniprot_mapping, query_pdb_code, query_chain_code)
    for query_uniprot in query_uniprots
        # AFDB
        filter!(m -> !occursin("AF-$(query_uniprot)-", m.target), search_results)
        # PDB
        query_structures = get_pdb_codes(sifts_uniprot_mapping, query_uniprot)
        for row in eachrow(query_structures)
            pdb = String(row.PDB)
            chain = String(row.CHAIN)
            filter!(m -> !_is_chain(String(m.target), pdb, chain), search_results)
        end
    end
    search_results
end

# Functions to look for known conformations of the proteins in the FoldSeek search results

# look for known conformations of the proteins showing similar structures to the query protein
function list_known_conformations(search_results::DataFrames.DataFrame, 
        sifts_uniprot_mapping::DataFrames.DataFrame)
    uniprots = Set{String}()
    new_targets = Set{String}()
    for row in eachrow(search_results)
        afdb_up = match(r"AF-(\w+)-F", row.target) # AFDB identifiers
        if afdb_up !== nothing
            push!(uniprots, afdb_up.captures[1])
        else # assume PDB if it is not AFDB
            pdb_code, chain_code = _get_pdb_and_chain(row.target)
            ups = get_uniprot_acc(sifts_uniprot_mapping, pdb_code, chain_code)
            for up in ups
                push!(uniprots, up)
            end
        end
    end
    for up in uniprots
        pdbs = get_pdb_codes(sifts_uniprot_mapping, up)
        for row in eachrow(pdbs)
            pdb = String(row.PDB)
            chain = String(row.CHAIN)
            # If the pdb is not in the search results, add it to the new targets
            if !any(target -> _is_chain(String(target), pdb, chain), search_results.target)
                push!(new_targets, "$(pdb).pdb_$(chain)") # Using the FoldSeek target format
            end
        end
    end
    new_targets
end

function add_known_conformations!(search_results::DataFrames.DataFrame, 
    sifts_uniprot_mapping::DataFrames.DataFrame)
    new_targets = list_known_conformations(search_results, sifts_uniprot_mapping)
    to_add = DataFrames.DataFrame(target=collect(new_targets))
    DataFrames.append!(search_results, to_add, cols=:union)
end