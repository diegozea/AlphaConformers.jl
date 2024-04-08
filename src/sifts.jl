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
function get_uniprot_mapping(sifts_db::String=get(ENV, "SIFTS_DB", pwd()))
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
    pdb::String, chain::Union{String,Type{MIToS.PDB.All}}=MIToS.PDB.All)
    pdb_code = lowercase(pdb)
    if chain !== MIToS.PDB.All
        ups = data[(data.PDB.==pdb_code).&(data.CHAIN.==chain), :SP_PRIMARY]
    else
        ups = data[data.PDB.==pdb_code, :SP_PRIMARY]
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
    pdb_codes = data[data.SP_PRIMARY.==uni_acc, [:PDB, :CHAIN]]
    uppercase.(unique(pdb_codes))
end

# FoldSeek results
# ================

# Functions to filter out the query protein from the FoldSeek search results

"""
    _get_pdb_and_chain(pdb_file::AbstractString)

Extracts the PDB code and chain identifier from a given `pdb_file` name. It returns
a tuple where the first element is the PDB code as an uppercase string and the second
element is the chain identifier. If the chain identifier is not provided, it defaults
to `MIToS.PDB.All`. 

!!! info
    This function only keeps uppercase letters and numbers from the `pdb_file` name. Then, 
    it assumes that the first 4 characters are the PDB code and the rest is the chain.
"""
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

"""
    _is_chain(pdb_file::String, pdb_code::String, chain_code::String)

Determines if a given `pdb_file` matches the specified `pdb_code` and `chain_code`. It does
so by extracting the actual PDB code and chain from `pdb_file` and comparing them to the
expected `pdb_code` and `chain_code`. Returns `true` if they match, `false` otherwise. This
function handles cases where the chain identifier might be missing (using `MIToS.PDB.All`).
"""
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

"""
    delete_query_from_target!(search_results::DataFrames.DataFrame,
                              sifts_uniprot_mapping::DataFrames.DataFrame,
                              query_pdb_code::String, query_chain_code::String)

Modifies `search_results` in place by removing any rows where the target matches the
query protein, identified by `query_pdb_code` and `query_chain_code`. This function
uses UniProt accession to identify all PDB codes and chains associated with the query
protein and filters out any corresponding matches from `search_results`. Both AFDB
and PDB formats are considered for matches.

!!! info
    This function is useful when testing the pipeline to avoid having the known 
    conformations of the query protein in the search results.
"""
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

"""
    uniprots_from_results(search_results::DataFrames.DataFrame,
                          sifts_uniprot_mapping::DataFrames.DataFrame)

It returns a Dict from the protein targets' identifiers on the `search_results` DataFrame
to UniProt accession numbers. It supports both AFDB and PDB formats for the identifiers.
For PDB entries, it translates PDB codes and chains to UniProt accessions using the
`sifts_uniprot_mapping` DataFrame—obtainable from SIFTS DB using the `get_uniprot_mapping`
function.
"""
function uniprots_from_results(search_results::DataFrames.DataFrame,
    sifts_uniprot_mapping::DataFrames.DataFrame)
    target2uniprot = Dict{String, String}()
    for row in eachrow(search_results)
        target = row.target
        afdb_up = match(r"AF-(\w+)-F", target) # AFDB identifiers
        if afdb_up !== nothing
            push!(target2uniprot, target => afdb_up.captures[1])
        else # assume PDB if it is not AFDB
            pdb_code, chain_code = _get_pdb_and_chain(target)
            ups = get_uniprot_acc(sifts_uniprot_mapping, pdb_code, chain_code)
            for up in ups
                push!(target2uniprot,  target => up)
            end
        end
    end
    target2uniprot
end

"""
    known_uniprot_structures!(target2uniprot::Dict{String,String},
                              sifts_uniprot_mapping::DataFrames.DataFrame,
                              search_results::DataFrames.DataFrame)::Set{String}

Identifies PDB codes and chains associated with a given set of UniProt accession numbers
(the unique values in `target2uniprot`) that are not already included in `search_results`. 
This function checks for new known structures by consulting the 
`sifts_uniprot_mapping` DataFrame. Returns a set of new PDB targets 
(formatted as "PDB.pdb_CHAIN") for inclusion in the search results. It also updates the
`target2uniprot` dictionary with the new targets.
"""
function known_uniprot_structures!(target2uniprot::Dict{String,String},
    sifts_uniprot_mapping::DataFrames.DataFrame,
    search_results::DataFrames.DataFrame)::Set{String}
    new_targets = Set{String}()
    uniprots = unique(values(target2uniprot))
    for up in uniprots
        pdbs = get_pdb_codes(sifts_uniprot_mapping, up)
        for row in eachrow(pdbs)
            pdb = String(row.PDB)
            chain = String(row.CHAIN)
            if !any(target -> _is_chain(String(target), pdb, chain), search_results.target)
                target_name = "$(pdb).pdb_$(chain)" # Using the FoldSeek target format
                push!(new_targets, target_name)
                push!(target2uniprot, target_name => up)
            end
        end
    end
    new_targets
end


# look for known conformations of the proteins showing similar structures to the query protein
"""
    list_known_conformations(search_results::DataFrames.DataFrame,
                             sifts_uniprot_mapping::DataFrames.DataFrame)

Compiles a list of known conformations of the proteins listed in the `search_results`. It
returns a tuple with a Dict from the protein targets' identifiers to UniProt accession and
a set of PDB codes and chains that are not currently included in the results.

!!! info
    This function first extracts UniProt accession numbers from `search_results`, 
    then identifies new PDB codes and chains associated with these accessions, 
    and returns these new targets.
"""
function list_known_conformations(search_results::DataFrames.DataFrame,
    sifts_uniprot_mapping::DataFrames.DataFrame)
    target2uniprot = uniprots_from_results(search_results, sifts_uniprot_mapping)
    new_targets = known_uniprot_structures!(
        target2uniprot, sifts_uniprot_mapping,  search_results)
    target2uniprot, new_targets
end

"""
    add_known_conformations!(search_results::DataFrames.DataFrame,
                             sifts_uniprot_mapping::DataFrames.DataFrame)

Expands `search_results` in place by adding new protein conformations identified as
related to the proteins in `search_results`. It utilizes `list_known_conformations`
to find these new conformations and appends them to the `search_results`. It return 
a tuple with a Dict from the protein targets' identifiers to UniProt accession and the
updated `search_results` table.
"""
function add_known_conformations!(search_results::DataFrames.DataFrame,
    sifts_uniprot_mapping::DataFrames.DataFrame)
    target2uniprot, new_targets = list_known_conformations(search_results, sifts_uniprot_mapping)
    to_add = DataFrames.DataFrame(target=collect(new_targets))
    DataFrames.append!(search_results, to_add, cols=:union)
    target2uniprot, search_results
end


"""
    get_uniprot2targets(target2uniprot::Dict{String,String}, 
                        search_results::DataFrames.DataFrame)

Given a `target2uniprot` dictionary mapping target structures to UniProt accessions,
this function returns a dictionary mapping UniProt accessions to a list of target 
structures. Note that only target structures from the original search results are
considered.
"""
function get_uniprot2targets(target2uniprot::Dict{String,String}, 
    search_results::DataFrames.DataFrame)
    original_targets = Set{String}(row.target for 
        row in eachrow(search_results) if !ismissing(row.evalue))
    uniprot2targets = Dict{String, Vector{String}}()
    for (target, uniprot) in target2uniprot
        if target in original_targets
            if haskey(uniprot2targets, uniprot)
                push!(uniprot2targets[uniprot], target)
            else
                uniprot2targets[uniprot] = [target]
            end
        end
    end
    uniprot2targets
end