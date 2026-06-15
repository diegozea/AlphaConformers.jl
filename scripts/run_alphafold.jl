# diego.zea@node48:/alpha/runs/diego.zea$ cp -r from_apo /store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers

const WORKING_DIR = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/from_apo"

const folders = ["1AEL_A", "1BV2_A", "1C54_A", "1EAL_A"]

function run_alphafold(folder)
    cd(WORKING_DIR)
    cd(folder)
    cd("clusters")
    subfolders = filter(f -> isdir(f) && startswith(f, "cluster_"), readdir())
    for subfolder in subfolders
        cd(subfolder)
        run(
            `/opt/alphafold/runcolabfold.py sequences.a3m af --use-templates 1 --msa-input --custom-template-path templates/ --num-seeds 5 --use-dropout --num-models 2`,
        )
        cd("..")
    end
    cd("..")
end

for folder in folders
    run_alphafold(folder)
end
