using AlphaConformers

const PATH = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset"
cd(PATH)

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

const FOLDSEEK_DB = "/alpha/database/pdb/fullpdb"

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"

const FOLDERS = filter!(isdir, readdir(join = true))

for folder in FOLDERS
    try
        run_alphafold(joinpath(folder, "clusters"), colabfold_path = COLABFOLD_PATH)
    catch err
        @error "Error processing $(folder): $(err)"
    end
end
