import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
import Foldseek_jll
import MIToS

function get_pdb(folder, pdb)
    pdb_file = joinpath(folder, "$(pdb).pdb")
    if !isfile(pdb_file)
        compressed_pdb_file = abspath("$(pdb_file).gz")
        MIToS.PDB.downloadpdb(pdb, filename=compressed_pdb_file, format=MIToS.PDB.PDBFile)
        res = MIToS.PDB.read(compressed_pdb_file, MIToS.PDB.PDBFile)
        rm(compressed_pdb_file)
        write(pdb_file, res, MIToS.PDB.PDBFile)
    end
end

# is the database folder already exists, delete all the database file
# while keeping the structures
function prepare_db_folder(folder)
    db_name = basename(folder)
    if isdir(folder)
        for file in readdir(folder)
            if startswith(file, db_name)
                rm(joinpath(folder, file))
            end
        end
    else
        mkdir(folder)
    end
end

const WORKING_DIR = @__DIR__
cd(WORKING_DIR)

const TEST_DB_PATH = joinpath(WORKING_DIR, "test_db")
prepare_db_folder(TEST_DB_PATH)

download_pdbs = false

if download_pdbs
    for pdb in ["4F4J", "1EX7", "3IIS"] # TRUE TRUE FALSE (query: 1EX6_B)
        get_pdb(TEST_DB_PATH, pdb)
    end
end

cd(TEST_DB_PATH)
run(`$(Foldseek_jll.foldseek()) createdb . test_db`)
cd(WORKING_DIR)

# Second database to check result merging
# ---------------------------------------

const TEST_DB_PATH_2 = joinpath(WORKING_DIR, "test_db_2")
prepare_db_folder(TEST_DB_PATH_2)

for pdb in ["1GKY", "8FXN"] # TRUE FALSE
    get_pdb(TEST_DB_PATH_2, pdb)
end

cd(TEST_DB_PATH_2)
run(`$(Foldseek_jll.foldseek()) createdb . test_db_2`)
cd(WORKING_DIR)