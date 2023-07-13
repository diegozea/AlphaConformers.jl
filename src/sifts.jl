# This file defines a series of utils to use the UniProt annotations from SIFT to 
# identify known conformations from a series of PDB files. Using the UniProt mapping 
# also allows to filter out structures from the query protein to evaluate the pipeline.

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
