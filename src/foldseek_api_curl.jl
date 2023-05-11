using MIToS.PDB
using JSON

pdb_code, chain_code = "1EX6", "B"

# pdb_file = downloadpdb(pdb_code, format=PDBFile)
chain = read("test/data/$(pdb_code).pdb.gz", PDBFile, chain=chain_code)
pdb_file = abspath("test/data/$(pdb_code)_$(chain_code).pdb")
write("test/data/$(pdb_code)_$(chain_code).pdb", chain, PDBFile)

function foldseek_api_query(pdb_file)
    foldseek_output = read(`curl -X POST -F q=@$pdb_file -F 'mode=3diaa' -F 'database[]=afdb50' -F 'database[]=afdb-swissprot' -F 'database[]=afdb-proteome' -F 'database[]=cath50' -F 'database[]=mgnify_esm30' -F 'database[]=pdb100' -F 'database[]=gmgcl_id' https://search.foldseek.com/api/ticket`, String)
    foldseek_output_json = JSON.parse(foldseek_output)
    ticket = foldseek_output_json["id"]
    if foldseek_output_json["status"] != "ERROR" && foldseek_output_json["status"] != "COMPLETE"
        status = JSON.parse(read(`curl -H "Content-Type: application/json" https://search.mmseqs.com/api/ticket/$ticket`, String))
        while status["status"] != "COMPLETE" && status["status"] != "ERROR"
            sleep(1)
            status = JSON.parse(read(`curl -H "Content-Type: application/json" https://search.mmseqs.com/api/ticket/$ticket`, String))
        end
        if status["status"] == "ERROR"
            throw(ErrorException("FoldSeek error"))
        end
    end
    read(`curl https://search.mmseqs.com/api/result/$ticket`, String)
end

foldseek_api_query(pdb_file)

# at the moment I am getting thsi error
#
# # MMseqs2 search will be back
# We are currently lacking resources to host both ColabFold and the MMseqs2 search server. 
# The search server will be back as soon as we can migrate ColabFold to different hardware.