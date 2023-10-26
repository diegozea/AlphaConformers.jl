import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
import Foldseek_jll
import MIToS

const WORKING_DIR = @__DIR__
cd(WORKING_DIR)

const TEST_DB_PATH = joinpath(WORKING_DIR, "test_db")
isdir(TEST_DB_PATH) || mkdir(TEST_DB_PATH)


download_pdbs = false

if download_pdbs
    for pdb in ["4F4J", "1EX7", "3IIS"] # TRUE TRUE FALSE (query: 1EX6_B)
        pdb_file = joinpath(TEST_DB_PATH, "$(pdb).pdb")
        if !isfile(pdb_file)
            compressed_pdb_file = joinpath("$(pdb_file).gz")
            MIToS.PDB.downloadpdb(pdb, filename=compressed_pdb_file, format=MIToS.PDB.PDBFile)
            res = MIToS.PDB.read(compressed_pdb_file, MIToS.PDB.PDBFile)
            rm(compressed_pdb_file)
            write(pdb_file, res, MIToS.PDB.PDBFile)
        end
    end
end

cd(TEST_DB_PATH)

run(`$(Foldseek_jll.foldseek()) createdb . test_db`)

cd(WORKING_DIR)
