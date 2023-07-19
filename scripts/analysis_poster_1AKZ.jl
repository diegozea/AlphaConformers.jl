using MIToS.PDB, AlphaConformers, CSV, DataFrames, Plots, StatsBase, Statistics

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

clusters_folder = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/1AKZ_A/clusters"
cd(clusters_folder)

cluster_folders = filter!(dir -> occursin("cluster_", dir) && "af" in readdir(dir), 
    readdir(join=true))

# 1AKZ_A;1SSP_E

apo_pdb = "1AKZ"
apo_chain = "A" 
holo_pdb = "1SSP"
holo_chain = "E"

# Figure 2. Apo (1AKZ chain A) and Holo (1SSP chain E) conformations of the Uracil-DNA glycosylase (UniProt ID: P13051)

apo_pdb_file = joinpath(PDB_FOLDER, "$apo_pdb.pdb")
holo_pdb_file = joinpath(PDB_FOLDER, "$holo_pdb.pdb")

apo_pdb_chain = read(apo_pdb_file, PDB.PDBFile, chain=apo_chain)
holo_pdb_chain = read(holo_pdb_file, PDB.PDBFile, chain=holo_chain)

apo_chain_file = "$(apo_pdb)_$(apo_chain).pdb"
holo_chain_file = "$(holo_pdb)_$(holo_chain).pdb"

write(apo_chain_file, apo_pdb_chain, PDB.PDBFile)
write(holo_chain_file, holo_pdb_chain, PDB.PDBFile)

apo_vs_holo =  usalign(apo_chain_file, holo_chain_file)

af_models = map(cluster_folders) do cl
    readdir(joinpath(cl, "af", "predictions", "sequences", "models"), join=true)
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

scatter(apo_vs_af_rmsd, holo_vs_af_rmsd, 
    xlabel="RMSD to apo", ylabel="RMSD to holo", 
    legend=false)
vline!([0.5295])
hline!([1.3836])

collect(Iterators.flatten(af_models))[findall(==(0), holo_vs_af_rmsd)]
holo = collect(Iterators.flatten(holo_vs_af))[findall(==(0), holo_vs_af_rmsd)]

holo_model = collect(Iterators.flatten(af_models))[findall(==(0.8), holo_vs_af_rmsd)]
holo_model_rmsd = collect(Iterators.flatten(holo_vs_af))[findall(==(0.8), holo_vs_af_rmsd)]

plot(read(first(holo_model), PDB.PDBFile))
plot!(read(holo_chain_file, PDB.PDBFile))

apo_model = collect(Iterators.flatten(af_models))[findall(==(0.5), apo_vs_af_rmsd)]
apo_model_rmsd = collect(Iterators.flatten(apo_vs_af))[findall(==(0.5), apo_vs_af_rmsd)]

cluster_numbers = map(af_models) do cl
    map(cl) do c
        @show c
        parse(Int, match(r"cluster_(\d+)", c).captures[1])
    end
end |> Iterators.flatten |> collect

sel = (holo_vs_af_rmsd .> 0) .& (apo_vs_af_rmsd .> 0)

gr(size=(600, 600))
max = maximum([maximum(apo_vs_af_rmsd), maximum(holo_vs_af_rmsd)])
ofs = 0.5
plt = scatter(apo_vs_af_rmsd[sel], holo_vs_af_rmsd[sel], 
    xlabel="RMSD to apo", ylabel="RMSD to holo", 
    marker_z=cluster_numbers[sel], markerstrokecolor=nothing, 
    aspect_ratio=1, label=nothing,
    markercolor=:viridis,
    xlims=(0, max + ofs), ylims=(0, max + ofs),
    ticks=0:1:(max + ofs))
scatter!(plt, [0.5295], [1.3836], c=:white, label=nothing)
hline!(plt, [1], c=:lightgrey, style=:dash, label=nothing)
vline!(plt, [1], c=:lightgrey, style=:dash, label=nothing)

Plots.svg(plt, "poster_1AKZ.svg")

# Figure 3. This figure shows the AlphaFold models (colored by cluster). The x-axis 
# shows the RMSD to the apo structure (1AKZ chain A) and the y-axis shows the RMSD to
# the holo structure (1SSP chain E). The white dot shows the AlphaFold model from 
# Saldaño et al. 2022 — default run of AlphaFold 2 using ColabFold. The PDB structures show 
# the supperposition of the AlphaFold models (blue) that are closer to the apo (oragne) and 
# holo (green) structures.