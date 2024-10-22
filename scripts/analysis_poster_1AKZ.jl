using MIToS.PDB, AlphaConformers, CSV, DataFrames, Plots, StatsBase, Statistics

const PDB_FOLDER = "/alpha/database/pdb/pdb_files"

# Apo	Holo
# 1AKZ_A	1SSP_E
# 1HW1_B	1H9G_A
# 1GQN_A	1QFE_B
# 2LHS-6_A	2BEM_A
# 1FMF-4_A	1ID8-11_A
# 2F63-4_A	1EQM_A
# 2UZ5-9_A	2VCD-8_A
# 4AKE_B	2ECK_B
# 6OY9_B	1A0R_B
# 1MUT-11_A	1PUN-7_A

# 1AKZ_A;1SSP_E # 0.5 0.8
# 1HW1_B;1H9G_A # 0.91 5.17
# 1GQN_A;1QFE_B # 0.47 0.47
# 2LHS-6_A;2BEM_A # 1.35 0.43
# 1FMF-4_A;1ID8-11_A # 1.7 3.01
# 2F63-4_A;1EQM_A # 2.09 1.55
# 2UZ5-9_A;2VCD-8_A # 2.14 1.57
# 4AKE_B;2ECK_B # 1.2 0.61
# 6OY9_B;1A0R_B # 1.37 1.49
# 1MUT-11_A;1PUN-7_A # 3.19 2.66

apo_dir = "1HW1_B"
apo = "1HW1_B"
holo = "1H9G_A"

#=
apo_dir = "2LHS-6_A"
apo = "2LHS_A"
holo = "2BEM_A"
=#

apo_pdb = String(split(apo, "_")[1])
apo_chain = String(split(apo, "_")[2])
holo_pdb = String(split(holo, "_")[1])
holo_chain = String(split(holo, "_")[2])

clusters_folder = "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/$(apo_dir)/clusters"
cd(clusters_folder)

cluster_folders = filter!(dir -> occursin("cluster_", dir) && "af" in readdir(dir), 
    readdir(join=true))


#=
apo_pdb = "1AKZ"
apo_chain = "A" 
holo_pdb = "1SSP"
holo_chain = "E"
=#

# Figure 2. Apo (1AKZ chain A) and Holo (1SSP chain E) conformations of the Uracil-DNA glycosylase (UniProt ID: P13051)

apo_pdb_file = joinpath(PDB_FOLDER, "$apo_pdb.pdb");
holo_pdb_file = joinpath(PDB_FOLDER, "$holo_pdb.pdb");

apo_pdb_chain = read(apo_pdb_file, PDB.PDBFile, chain=apo_chain);
holo_pdb_chain = read(holo_pdb_file, PDB.PDBFile, chain=holo_chain);

apo_chain_file = "$(apo_pdb)_$(apo_chain).pdb";
holo_chain_file = "$(holo_pdb)_$(holo_chain).pdb";

write(apo_chain_file, apo_pdb_chain, PDB.PDBFile);
write(holo_chain_file, holo_pdb_chain, PDB.PDBFile);

apo_vs_holo =  usalign(apo_chain_file, holo_chain_file);

af_models = map(cluster_folders) do cl
    path = joinpath(cl, "af", "predictions", "sequences", "models")
    if isdir(path)
        readdir(path, join=true)
    else
        String[]
    end
end;

apo_vs_af = map(af_models) do cluster
    map(cluster) do af
        usalign(apo_chain_file, af)
    end
end;

holo_vs_af = map(af_models) do cluster
    map(cluster) do af
        usalign(holo_chain_file, af)
    end
end;

apo_vs_af_rmsd = map(apo_vs_af) do cluster
    map(cluster) do af
        only(af.RMSD)
    end
end |> Iterators.flatten |> collect;

holo_vs_af_rmsd = map(holo_vs_af) do cluster
    map(cluster) do af
        only(af.RMSD)
    end
end |> Iterators.flatten |> collect;


minimum(filter(>(0), apo_vs_af_rmsd))
minimum(filter(>(0), holo_vs_af_rmsd))



#=
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

=#