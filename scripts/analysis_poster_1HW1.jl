using MIToS.PDB, AlphaConformers, CSV, DataFrames, Plots, StatsBase, Statistics

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

clusters_folder = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/1HW1_B/clusters"
cd(clusters_folder)

cluster_folders =
    filter!(dir -> occursin("cluster_", dir) && "af" in readdir(dir), readdir(join = true))

apo_pdb = "1HW1"
apo_chain = "B"
holo_pdb = "1H9G"
holo_chain = "A"

apo_pdb_file = joinpath(PDB_FOLDER, "$apo_pdb.pdb")
holo_pdb_file = joinpath(PDB_FOLDER, "$holo_pdb.pdb")

apo_pdb_chain = read(apo_pdb_file, PDB.PDBFile, chain = apo_chain)
holo_pdb_chain = read(holo_pdb_file, PDB.PDBFile, chain = holo_chain)

apo_chain_file = "$(apo_pdb)_$(apo_chain).pdb"
holo_chain_file = "$(holo_pdb)_$(holo_chain).pdb"

write(apo_chain_file, apo_pdb_chain, PDB.PDBFile)
write(holo_chain_file, holo_pdb_chain, PDB.PDBFile)

apo_vs_holo = usalign(apo_chain_file, holo_chain_file)

af_models = map(cluster_folders) do cl
    readdir(joinpath(cl, "af", "predictions", "sequences", "models"), join = true)
end

apo_vs_af = map(af_models) do cluster
    map(cluster) do af
        usalign(apo_chain_file, af)
    end
end

holo_vs_af = map(af_models) do cluster
    map(cluster) do af
        usalign(holo_chain_file, af)
    end
end

apo_vs_af_rmsd = map(apo_vs_af) do cluster
    map(cluster) do af
        only(af.RMSD)
    end
end |> Iterators.flatten |> collect

holo_vs_af_rmsd = map(holo_vs_af) do cluster
    map(cluster) do af
        only(af.RMSD)
    end
end |> Iterators.flatten |> collect

plotly()

scatter(
    apo_vs_af_rmsd,
    holo_vs_af_rmsd,
    xlabel = "RMSD to apo",
    ylabel = "RMSD to holo",
    legend = false,
)

collect(Iterators.flatten(af_models))[findall(==(0), holo_vs_af_rmsd)]
holo = collect(Iterators.flatten(holo_vs_af))[findall(==(0), holo_vs_af_rmsd)]
