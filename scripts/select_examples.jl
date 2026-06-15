using CSV, DataFrames, Plots
plotly()

data = DataFrame(
    CSV.File(
        "/stockage/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/AlphaConformers/data/Supplementary_Table_1_91_apo_holo_pairs.csv",
    ),
)

scatter(
    data.RMSD_best_score_md_Apo,
    data.RMSD_best_score_md_Holo,
    xlabel = "RMSD best score apo",
    ylabel = "RMSD best score holo",
    legend = false,
    hover = 1:nrow(data),
)
plot!(y -> y, hover = "")
hline!([1.0], hover = "")
vline!([1.0], hover = "")

# to select: 
# 56, 73 (good apo, acceptable holo)
# 8, 36 (good holo, acceptable apo)
# 6, 27 (acceptable apo, bad holo)
# 39, 40 (bad apo, acceptable holo)
# 19, 82 (both acceptable)
# 16, 88 (both bad)

to_select = [56, 73, 8, 36, 6, 27, 39, 40, 19, 82, 16, 88]

subset = data[to_select, :]

poster_subset_path = "/stockage/EQUIPES/AMIG/MEMBERS/diego.zea/AlphaConformers/poster_subset"
CSV.write(joinpath(poster_subset_path, "selected_examples.csv"), subset)
