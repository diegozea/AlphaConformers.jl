using Distributed
addprocs(12)

@everywhere using MIToS.PDB, AlphaConformers, CSV, DataFrames, Plots, StatsBase, Statistics
@everywhere import UnicodePlots
@everywhere import MIToS
@everywhere import JSON3

@everywhere const PDB_FOLDER = "/alpha/database/pdb/pdb_files"
!isdir(PDB_FOLDER) && @warn "PDB database is located at node 48"

@everywhere const EXAMPLES = DataFrame(
    CSV.File(
        "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/selected_examples.csv",
    ),
)

@everywhere function get_path_clusters_folder(apo_pdb, apo_chain)
    "/store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/$(apo_pdb)_$(apo_chain)/clusters"
end

@everywhere function clean_and_get_model_id(code::String)
    fields = split(code, "-")
    if length(fields) == 1
        return code, MIToS.PDB.All
    else
        return String(fields[1]), String(fields[2])
    end
end

@everywhere function compare_models(apo_pdb, apo_chain, holo_pdb, holo_chain)
    clusters_folder = get_path_clusters_folder(apo_pdb, apo_chain)
    cwd = pwd()
    cd(clusters_folder)
    apo_vs_holo = nothing
    af_models = nothing
    apo_vs_af = nothing
    holo_vs_af = nothing
    apo_vs_af_rmsd = nothing
    holo_vs_af_rmsd = nothing
    scores = nothing
    msas = nothing
    try
        cluster_folders = filter!(
            dir -> occursin("cluster_", dir) && "af" in readdir(dir),
            readdir(join = true),
        )

        apo_chain_file = "$(apo_pdb)_$(apo_chain).pdb"
        holo_chain_file = "$(holo_pdb)_$(holo_chain).pdb"

        if !isfile(apo_chain_file) || !isfile(holo_chain_file)
            clean_apo_pdb, apo_model = clean_and_get_model_id(apo_pdb)
            clean_holo_pdb, holo_model = clean_and_get_model_id(holo_pdb)

            apo_pdb_file = joinpath(PDB_FOLDER, "$clean_apo_pdb.pdb")
            holo_pdb_file = joinpath(PDB_FOLDER, "$clean_holo_pdb.pdb")

            apo_pdb_chain =
                read(apo_pdb_file, PDB.PDBFile, chain = apo_chain, model = apo_model)
            holo_pdb_chain =
                read(holo_pdb_file, PDB.PDBFile, chain = holo_chain, model = holo_model)

            write(apo_chain_file, apo_pdb_chain, PDB.PDBFile)
            write(holo_chain_file, holo_pdb_chain, PDB.PDBFile)
        end

        apo_vs_holo = usalign(apo_chain_file, holo_chain_file)

        af_models = []
        msas = MIToS.MSA.AnnotatedMultipleSequenceAlignment[]
        for cl in cluster_folders
            path = joinpath(cl, "af", "predictions", "sequences", "models")
            if isdir(path)
                push!(af_models, readdir(path, join = true))
                push!(
                    msas,
                    MIToS.MSA.read(joinpath(cl, "af", "sequences.a3m"), MIToS.MSA.FASTA),
                )
            end
        end

        scores = map(af_models) do cluster
            map(cluster) do af
                scores_path = replace(
                    replace(replace(af, ".pdb" => ".json"), "_unrelaxed_" => "_scores_"),
                    "sequences/models" => "sequences/scores",
                )
                JSON3.read(read(scores_path, String))
            end
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
    catch err
        rethrow(err)
    finally
        cd(cwd)
    end
    (;
        apo_vs_holo,
        af_models,
        apo_vs_af,
        holo_vs_af,
        apo_vs_af_rmsd,
        holo_vs_af_rmsd,
        scores,
        msas,
    )
end

@everywhere function plot_data(comparison)
    holo_vs_af_rmsd = comparison.holo_vs_af_rmsd
    apo_vs_af_rmsd = comparison.apo_vs_af_rmsd

    sel = (holo_vs_af_rmsd .> 0) .& (apo_vs_af_rmsd .> 0)

    plt = UnicodePlots.scatterplot(
        apo_vs_af_rmsd[sel],
        holo_vs_af_rmsd[sel],
        xlabel = "RMSD to apo",
        ylabel = "RMSD to holo",
    )
    UnicodePlots.hline!(plt, 1.5)
    UnicodePlots.vline!(plt, 1.5)
    plt
end

# julia> EXAMPLES[:, 1:2]
# 12×2 DataFrame
#  Row │ apo_id     holo_id   
#      │ String15   String15  
# ─────┼──────────────────────
#    1 │ 1AKZ_A     1SSP_E
#    2 │ 1HW1_B     1H9G_A
#    3 │ 1GQN_A     1QFE_B
#    4 │ 2LHS-6_A   2BEM_A
#    5 │ 1FMF-4_A   1ID8-11_A
#    6 │ 2F63-4_A   1EQM_A
#    7 │ 2UZ5-9_A   2VCD-8_A
#    8 │ 4AKE_B     2ECK_B
#    9 │ 1O1U-9_A   1O1V-5_A
#   10 │ 6OY9_B     1A0R_B
#   11 │ 1MUT-11_A  1PUN-7_A
#   12 │ 4LP5_A     4P2Y_A

# 56, 73 (good apo, acceptable holo)
# 8, 36 (good holo, acceptable apo)
# 6, 27 (acceptable apo, bad holo)
# 39, 40 (bad apo, acceptable holo)
# 19, 82 (both acceptable)
# 16, 88 (both bad)

# comparison = compare_models("1AKZ", "A", "1SSP", "E");      ## APO OK HOLO OK : + + +  # (good apo, acceptable holo) # domain movement
# comparison = compare_models("1HW1", "B", "1H9G", "A");      ## APO OK HOLO NO : - + -  # (good apo, acceptable holo) # domain movement
# comparison = compare_models("1GQN", "A", "1QFE", "B");      ## APO OK HOLO OK : + + +  # (good holo, acceptable apo) # loop movement
# comparison = compare_models("2LHS-6", "A", "2BEM", "A");    ## APO NO HOLO OK : - + +  # (good holo, acceptable apo) # loop movement
# comparison = compare_models("1FMF-4", "A", "1ID8-11", "A"); ## APO NO HOLO NO : - - -  # (acceptable apo, bad holo)  # loop movement
# comparison = compare_models("2F63-4", "A", "1EQM", "A");    ## APO NO HOLO NO : - - -  # (acceptable apo, bad holo)  # loop movement
# comparison = compare_models("2UZ5-9", "A", "2VCD-8", "A");  ## APO NO HOLO NO : - - -  # (bad apo, acceptable holo)  # loop movement
# comparison = compare_models("4AKE", "B", "2ECK", "B");      ## APO NO HOLO OK : - + +  # (bad apo, acceptable holo)  # domain movement
# comparison = compare_models("1O1U-9", "A", "1O1V-5", "A");  ## NO DATA                # (both acceptable)           # loop movement
# comparison = compare_models("6OY9", "B", "1A0R", "B");      ## APO NO HOLO NO : - + -  # (both acceptable)           # loop movement
# comparison = compare_models("1MUT-11", "A", "1PUN-7", "A"); ## APO NO HOLO NO : - - -  # (both bad)                  # loop movement
# comparison = compare_models("4LP5", "A", "4P2Y", "A");      ## ERROR                  # (both bad)                  # domain movement

# RERUN = ["1FMF-4_A", "2F63-4_A", "1O1U-9_A", "4LP5_A"] 

comparison = compare_models("1AKZ", "A", "1SSP", "E");
plot_data(comparison)

# --- all --- #

comparisons = pmap(eachrow(EXAMPLES)) do row
    apo_id = row.apo_id
    holo_id = row.holo_id
    apo_fields = split(apo_id, "_")
    holo_fields = split(holo_id, "_")
    apo_pdb, apo_chain = apo_fields[1], apo_fields[2]
    holo_pdb, holo_chain = holo_fields[1], holo_fields[2]
    try
        compare_models(apo_pdb, apo_chain, holo_pdb, holo_chain)
    catch err
        println("ERROR: $(apo_id) $(holo_id): $err")
        missing
    end
end;

# --- at frontend --- #

# comparison = compare_models("4AKE", "B", "2ECK", "B");
# comparison = compare_models("1AKZ", "A", "1SSP", "E");
# comparison = compare_models("2UZ5-9", "A", "2VCD-8", "A");
# comparison = compare_models("2LHS-6", "A", "2BEM", "A");
# comparison = compare_models("6OY9", "B", "1A0R", "B");
# "1FMF-4_A", "2F63-4_A"
# comparison = compare_models("1FMF-4", "A", "1ID8-11", "A");
# comparison = compare_models("2F63-4", "A", "1EQM", "A");

ptm = [[s.ptm for s in l] for l in comparison.scores] |> Iterators.flatten |> collect

plotly()

plot(
    scatter(comparison.apo_vs_af_rmsd, ptm, xlabel = "RMSD to apo", ylabel = "pTM"),
    scatter(comparison.holo_vs_af_rmsd, ptm, xlabel = "RMSD to holo", ylabel = "pTM"),
    legend = false,
)

comparison.msas

nseqs = MIToS.MSA.nsequences.(comparison.msas)
mean_pid = MIToS.MSA.meanpercentidentity.(comparison.msas)

gr()

plot(
    scatter(comparison.apo_vs_af_rmsd, nseqs, xlabel = "RMSD to apo", ylabel = "nseqs"),
    scatter(comparison.holo_vs_af_rmsd, nseqs, xlabel = "RMSD to holo", ylabel = "nseqs"),
    legend = false,
)

plot(
    scatter(
        comparison.apo_vs_af_rmsd,
        mean_pid,
        xlabel = "RMSD to apo",
        ylabel = "mean pid",
    ),
    scatter(
        comparison.holo_vs_af_rmsd,
        mean_pid,
        xlabel = "RMSD to holo",
        ylabel = "mean pid",
    ),
    legend = false,
)

plot(
    scatter(ptm, nseqs, xlabel = "pTM", ylabel = "nseqs"),
    scatter(ptm, mean_pid, xlabel = "pTM", ylabel = "mean pid"),
    legend = false,
)

# --- look for errors --- #

clusters_path = get_path_clusters_folder("6OY9", "B")

# ---
# 2LHS-6_A has sequences in the but there is no pdb in the template folder
# ; cat /store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/2LHS-6_A/clusters/cluster_2/sequences.fasta
#
# # # lowercase chain n in 7FJE.pdb
# # ┌ Warning: The are no ATOM residues in /store/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset/tmp/2LHS-6_A/pdb/7FJE.pdb_N (model: 1, chain: N)
# # └ @ AlphaConformers ~/.julia/dev/AlphaConformers/src/af_input.jl:69
#
# ---
# 60Y9_B has only the sequence file for a single cluster
# ---


# run(`$(Foldseek_jll.foldseek()) easy-search 2HBB.pdb /alpha/database/afdb/afdb_up output tmp`);
# Try: 35Gb RAM -> --prefilter-mode 1
# run(`$(Foldseek_jll.foldseek()) easy-search 2HBB.pdb /alpha/database/afdb/afdb_up output tmp --prefilter-mode 1`);
# NOTE: Using --prefilter-mode 1, Foldseek requires 35Gb RAM instead of 325Gb :)
# run(`$(Foldseek_jll.foldseek()) easy-search 2HBB.pdb /alpha/database/afdb/afdb_up output2 tmp2 --prefilter-mode 1 --format-mode 5`);


# --- at frontend --- #
# Question 1 : Are the generated models close to their templates?

# To answer this question, we will compare the generated models from a given cluster to all 
# the templates structures. That will allow us to visualize if the generated models are 
# closer to its cluster templates than to the other templates (that last case could mean
# that AlphaFold is ignoring the templates).

# Question 2 : Are problem arrising from using the whole protein chain rather than 
# delimiting to the aligned region?

# To answer this question, we will check wether the difference in legth of query and target
# proteins are correlated with the RMSD of the generated models. As I have the PDB files
# I will use the length of the PDB chain as a proxy for the protein length.

# 4AKE : The apo structure can not be modelled with AlphaFold

apo_pdb_code = "4AKE"
apo_chain_code = "B"
holo_pdb_code = "2ECK"
holo_chain_code = "B"
comparison = compare_models(apo_pdb_code, apo_chain_code, holo_pdb_code, holo_chain_code);
clusters_path = get_path_clusters_folder(apo_pdb_code, apo_chain_code)

readdir(joinpath(clusters_path, "..", "pdb"))
readdir(joinpath(clusters_path, "cluster_1", "templates"))

# It looks like AlphaFold is using a sequence alignment to design the templates anyway; so 
# it should be fine with AlphaFold. However, it is not clear if the MAFFT alignments are 
# correct in such cases. Maybe, I an also look at the gap fraction in the alignment for 
# example to have a proxy about its quality.

gaps = MIToS.MSA.gapfraction.(comparison.msas)

scatter(ptm, gaps, xlabel = "pTM", ylabel = "gap fraction")
scatter(comparison.apo_vs_af_rmsd, gaps, xlabel = "RMSD to apo", ylabel = "gap fraction")
scatter(comparison.holo_vs_af_rmsd, gaps, xlabel = "RMSD to holo", ylabel = "gap fraction")

# 4AKE : There are acceptable models for the apo structure with low and high gap fractions.

# 1: align templates to the apo structure
cd(clusters_path)
apo_pdb = "$(apo_pdb_code)_$(apo_chain_code).pdb"
@assert isfile(apo_pdb)
holo_pdb = "$(holo_pdb_code)_$(holo_chain_code).pdb"
@assert isfile(holo_pdb)

function cluster_stats(pdb, clusters_path, cluster_number)
    na = (
        min_rmsd = missing,
        mean_rmsd = missing,
        max_rmsd = missing,
        mean_len_diff = missing,
    )
    cluster_folder = joinpath(clusters_path, "cluster_$(cluster_number)", "templates")
    !isdir(cluster_folder) && return na
    template_files = filter!(endswith(".pdb"), readdir(cluster_folder))
    isempty(template_files) && return na
    aln_data = usalign(pdb, cluster_folder, template_files)
    aligned = aln_data[aln_data.Lali .> 20, :]
    len_diff = abs.(aligned.L2 .- aligned.L1)
    min_rmsd = minimum(aligned.RMSD)
    mean_rmsd = mean(aligned.RMSD)
    max_rmsd = maximum(aligned.RMSD)
    mean_len_diff = mean(len_diff)
    (; min_rmsd, mean_rmsd, max_rmsd, mean_len_diff)
end

apo_stats = cluster_stats(apo_pdb, clusters_path, 1)
holo_stats = cluster_stats(holo_pdb, clusters_path, 1)

cluster_numbers =
    parse.(
        Int,
        last.(
            split.(filter!(f -> isdir(f) && occursin('_', f), readdir(clusters_path)), '_'),
        ),
    ) |> sort

stats = map(cluster_numbers) do cl
    apo_stats = cluster_stats(apo_pdb, clusters_path, cl)
    holo_stats = cluster_stats(holo_pdb, clusters_path, cl)
    (; apo = apo_stats, holo = holo_stats)
end

min_apo = [
    s.apo.min_rmsd for
    s in stats if !ismissing(s.apo.min_rmsd) && !ismissing(s.holo.min_rmsd)
]
min_holo = [
    s.holo.min_rmsd for
    s in stats if !ismissing(s.apo.min_rmsd) && !ismissing(s.holo.min_rmsd)
]
scatter(min_apo, min_holo, xlabel = "min apo", ylabel = "min holo", legend = false)

# Well, in fact, there are not good templates for the apo structure.

# From comparison get the minimum RMSD for each cluster by looking into af_models and 
# apo_vs_af_rmsd/holo_vs_af_rmsd

data = DataFrame(
    cluster = [
        parse(Int, split(p[11], "_")[end]) for
        p in comparison.af_models |> Iterators.flatten |> collect .|> splitpath
    ],
    rmsd = comparison.apo_vs_af_rmsd,
)

data = data[data.rmsd .> 0, :]

sorted = combine(groupby(data, :cluster), :rmsd => minimum) |> sort

stats_data = DataFrame(
    cluster = cluster_numbers,
    min_apo = [s.apo.min_rmsd for s in stats],
    min_holo = [s.holo.min_rmsd for s in stats],
)

stats_data = stats_data[completecases(stats_data), :]

combined = innerjoin(stats_data, sorted, on = :cluster)

scatter(
    combined.min_apo,
    combined.min_holo,
    xlabel = "min apo",
    ylabel = "min holo",
    legend = false,
)

plotly()
plt = scatter(
    combined.min_apo,
    combined.rmsd_minimum,
    xlabel = "min apo templates",
    ylabel = "min apo models",
    legend = false,
)
plot!(plt, x -> x)

# Not having good templates for the apo structure looks to be the limiting factor for
# AlphaFold to model the apo structure. However, we arrive to have 1.2 RMSD for the apo
# in a cluster were the minimum RMSD for the templates to the apo structure is 1.5.

# Now; let's try to see why we do not get a good model for the holo structure in the 
# 1HW1_B - 1H9G_A pair.
apo_pdb_code = "1HW1"
apo_chain_code = "B"
holo_pdb_code = "1H9G"
holo_chain_code = "A"
comparison = compare_models(apo_pdb_code, apo_chain_code, holo_pdb_code, holo_chain_code);
clusters_path = get_path_clusters_folder(apo_pdb_code, apo_chain_code)

gaps = MIToS.MSA.gapfraction.(comparison.msas)

gr()
scatter(ptm, gaps, xlabel = "pTM", ylabel = "gap fraction")
scatter(comparison.apo_vs_af_rmsd, gaps, xlabel = "RMSD to apo", ylabel = "gap fraction")
scatter(comparison.holo_vs_af_rmsd, gaps, xlabel = "RMSD to holo", ylabel = "gap fraction")

cd(clusters_path)
apo_pdb = "$(apo_pdb_code)_$(apo_chain_code).pdb"
@assert isfile(apo_pdb)
holo_pdb = "$(holo_pdb_code)_$(holo_chain_code).pdb"
@assert isfile(holo_pdb)


cluster_numbers =
    parse.(
        Int,
        last.(
            split.(filter!(f -> isdir(f) && occursin('_', f), readdir(clusters_path)), '_'),
        ),
    ) |> sort

stats = map(cluster_numbers) do cl
    apo_stats = cluster_stats(apo_pdb, clusters_path, cl)
    holo_stats = cluster_stats(holo_pdb, clusters_path, cl)
    (; apo = apo_stats, holo = holo_stats)
end

min_apo = [
    s.apo.min_rmsd for
    s in stats if !ismissing(s.apo.min_rmsd) && !ismissing(s.holo.min_rmsd)
]
min_holo = [
    s.holo.min_rmsd for
    s in stats if !ismissing(s.apo.min_rmsd) && !ismissing(s.holo.min_rmsd)
]

p = scatter(min_apo, min_holo, xlabel = "min apo", ylabel = "min holo", legend = false)
plot!(p, x -> x)


data = DataFrame(
    cluster = [
        parse(Int, split(p[11], "_")[end]) for
        p in comparison.af_models |> Iterators.flatten |> collect .|> splitpath
    ],
    rmsd = comparison.apo_vs_af_rmsd,
)

data = data[data.rmsd .> 0, :]

sorted = combine(groupby(data, :cluster), :rmsd => minimum) |> sort

stats_data = DataFrame(
    cluster = cluster_numbers,
    min_apo = [s.apo.min_rmsd for s in stats],
    min_holo = [s.holo.min_rmsd for s in stats],
)

stats_data = stats_data[completecases(stats_data), :]

combined = innerjoin(stats_data, sorted, on = :cluster)

pp = scatter(
    combined.min_apo,
    combined.min_holo,
    xlabel = "min apo",
    ylabel = "min holo",
    legend = false,
)
plot!(pp, x -> x)

# plotly()
plt = scatter(
    combined.min_apo,
    combined.rmsd_minimum,
    xlabel = "min apo templates",
    ylabel = "min apo models",
    legend = false,
)
plot!(plt, x -> x)

plot(pp, plt)

# HOLO
data = DataFrame(
    cluster = [
        parse(Int, split(p[11], "_")[end]) for
        p in comparison.af_models |> Iterators.flatten |> collect .|> splitpath
    ],
    rmsd = comparison.holo_vs_af_rmsd,
)
data = data[data.rmsd .> 0, :]
sorted = combine(groupby(data, :cluster), :rmsd => minimum) |> sort
stats_data = DataFrame(
    cluster = cluster_numbers,
    min_apo = [s.apo.min_rmsd for s in stats],
    min_holo = [s.holo.min_rmsd for s in stats],
)
stats_data = stats_data[completecases(stats_data), :]
combined = innerjoin(stats_data, sorted, on = :cluster)
plt_holo = scatter(
    combined.min_holo,
    combined.rmsd_minimum,
    xlabel = "min holo templates",
    ylabel = "min holo models",
    legend = false,
)
plot!(plt_holo, x -> x)

plot(pp, plt, plt_holo, layout = (1, 3), link = :both)

# (apo = (min_rmsd = 0.42, mean_rmsd = 4.396833333333333, max_rmsd = 6.44, mean_len_diff = 83.06666666666666), holo = (min_rmsd = 0.39, mean_rmsd = 4.24214953271028, max_rmsd = 6.12, mean_len_diff = 88.05607476635514))

# A cluster (cluster 2) with a template close to the apo and holo structures, gives a model 
# that is far from both structures. It could be a proble with legth or with the structural 
# variability in the cluster as the mean RMSD is high between the structures and the 
# templates. I think atht maybe, since I am not checking the coverage when clutering; I am 
# getting good RMSD (lower than 1) but for a very small region. Therefore, this bus is 
# poluting the clusters.
