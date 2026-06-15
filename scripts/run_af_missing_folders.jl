using Revise
using AlphaConformers

const PATH = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset"
cd(PATH)

const COLABFOLD_PATH = "/opt/alphafold/runcolabfold.py"

const FOLDERS = abspath.([
    "1FMF-4_A",
    "1GQN_A",
    "1MUT-11_A",
    "1O1U-9_A",
    "2F63-4_A",
    "2LHS-6_A",
    "2UZ5-9_A",
    "4AKE_B",
    "4LP5_A",
    "6OY9_B",
])

for folder in FOLDERS[3:end]
    try
        run_alphafold(joinpath(folder, "clusters"), colabfold_path = COLABFOLD_PATH)
    catch err
        @warn "Error processing $(folder): $(err)"
    end
end
