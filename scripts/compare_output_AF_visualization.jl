import Pkg
Pkg.activate("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/update")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/data/") 
# Load necessary packages on all workers

using CSV
using DataFrames
using StatsPlots
using Plots
using Measures
using Statistics

gr()  # ou pyplot() selon ce que tu utilises

########################################################## FUNCTIONS #################################################
"""
Create the plot to compare the performance of ALphaConformer depending on database, evalue and cutoff
Have one plot for each pdb with all the performance for each parameter
Take in input the DF output by automatic_AlphaConformer_Test_results.jl
"""
#function to compare the parameters of AlphaConformer
function compare_parameters(df::DataFrame)

    @show df
    # get the unique names
    pdbs = unique(df.Query)

    # for each pdb, create a figure with the different parameters
    for pdb in pdbs
        subdf = filter(row -> row.Query == pdb, df)
        # get the unique values for database, evalue and cutoff
        databases = unique(subdf.Database)
        evalues = unique(subdf.Evalue)
        cutoffs = unique(subdf.Cutoff)
        
        # prepare the plot
        nrows = length(databases)
        ncols = length(evalues)
        plt = plot(
            layout = (nrows, ncols),
            size = (1400, 800),
            suptitle = "Query : $pdb with Foldseek", 
            bottom_margin = 10mm,
            top_margin = 5mm,
            right_margin = 5mm,
            left_margin = 5mm,
            legend = :outerright 
            )

        for (i, db) in enumerate(databases)
            for (j, ev) in enumerate(evalues)
                # Filteer the DataFrame for the current database and evalue
                sdf = filter(row -> row.Database == db && ((isnan(row.Evalue) && isnan(ev)) || row.Evalue == ev), subdf)
                # get the min_apo and min_holo values
                min_apo_vals = sdf.min_apo
                min_holo_vals = sdf.min_holo
                rmsd_val = first(sdf.rmsd_apo_holo)
                all_vals=[]
                for i in 1:length(min_apo_vals)
                    push!(all_vals,min_apo_vals[i])
                    push!(all_vals,min_holo_vals[i])
                end
                #data for the bar plot
                subplot_idx = (i - 1) * ncols + j
                cutoffs= unique(sdf.Cutoff)
                xticks_pos = [2i-0.5 for i in 1:length(cutoffs)]
                xticks_labels = string.(cutoffs)
                x_apo = 1:2:length(all_vals)
                x_holo = 2:2:length(all_vals)
                all_rmsd_vals = [first(filter(row -> row.Database == db && ((isnan(row.Evalue) && isnan(ev)) || row.Evalue == ev), subdf).rmsd_apo_holo)
                    for db in databases, ev in evalues]
                max_rmsd_val = maximum(all_rmsd_vals)
                ylim_global = (0, 2 * max_rmsd_val)

                # Create the bar plot
                bar!(
                    plt,
                    x_apo,
                    min_apo_vals;
                    bar_width = 0.8,
                    color = :orange,
                    label = "min apo",
                    subplot = subplot_idx,
                    legend = subplot_idx == 1,
                    title = "DB: $db, Evalue: $ev",
                    xlabel= "Cutoff",
                    ylabel = "RMSD",
                    ylim = ylim_global,
                    xticks = (xticks_pos, xticks_labels) 
                )
                bar!(
                    plt,
                    x_holo,
                    min_holo_vals;
                    bar_width = 0.8,
                    color = :lightgreen,
                    label = "min holo",
                    subplot = subplot_idx,
                    legend = subplot_idx == 1,
                    ylim = ylim_global
                )
                hline!(
                    plt,
                    [rmsd_val];
                    color = :red,
                    lw = 2,
                    linestyle = :dash,
                    label = "RMSD between apo et holo",
                    subplot = subplot_idx,
                    legend = subplot_idx == 1
                )

            end
        end
        # Save the plot
        savefig(plt, "comparaison_resultat_AlphaConformer_$pdb.png")
    end
end

"""
Visualize the min apo and min holo for each methods 
Put a red line for the RMSD between apo and holo
Take in input the DF output by compare_output_AF.jl
"""
# Compare the result between different pipeline 
function compare_methods(df)

    # Get the unique pdb names
    pdb_names = ["6akm","6uui","7lp1","9g21"] # can to it for all the file 
    #pdb_names = unique(df.Query)

    # Difine the colors for the bar plots
    colors_list = [
        :orange,          # 1er élément
        :lightgreen,    # 2
        :orange,          # 1er élément
        :lightgreen, 
        :orange,          # 1er élément
        :lightgreen,         # 7
    ]
    line_color = :red

    #for each pdb, create a figure with the different methods
    for pdb in pdb_names
        # Filter the DataFrame for the current pdb
        df_subset = filter(row -> row.File == pdb &&
        (row.Folder == pdb || row.Folder == pdb * "_AF_CLUSTER" || row.Folder == pdb * "_BIOEMU"), df)
        n = nrow(df_subset)
        labels = String.(df_subset.File)
        #get the min_apo and min_holo values
        min_apo_vals = df_subset.min_apo
        min_holo_vals = df_subset.min_holo
        rmsd_val = first(df_subset.rmsd_apo_holo)
        all_vals=[]
        println(min_apo_vals)
        for i in 1:length(min_apo_vals)
            push!(all_vals,min_apo_vals[i])
            push!(all_vals,min_holo_vals[i])
        end
        println(all_vals)
        # Barplot simple
        xticks_pos = [2i-0.5 for i in 1:length(all_vals)]
        xticks_labels = ["AlphaConformer", "AF-Cluster", "BioEmu"]
        p = bar(
            1:(2n),
            all_vals,
            bar_width = 0.8,
            c = colors_list[1:2n],  # Trim la liste à la bonne taille
            legend = false,
            xticks = (xticks_pos, xticks_labels),
            ylabel = "RMSD",
            title = "Comparaison RMSD pour $pdb",
            framestyle = :box,
            titlefont = font(14, "Arial"),
            tickfont = font(10)
        )

        
        # Ligne horizontale
        hline!(p, [rmsd_val], color = line_color, lw = 2, linestyle = :dash, label = "rmsd_apo_holo")

        savefig(p, "comparaison_methodes_$pdb.png")
    end

end

"""
Visualize the min apo and min holo for each methods and N pdbs
Take in input the DF output by compare_output_AF.jl
Can specify the pdb to compare 
Have one plot to compare all the pdb with the 3 pipeline 
"""
#Compare the result between different pipeline for multiple pdb in one plot
function compare_method_2(df)
    
    pdb_names = ["6akm","6uui","7lp1","9g21"] # can to it for all the file 
    #pdb_names = unique(df.Query)

    method_labels = String[]
    value_types = String[]
    values = Float64[]
    for pdb in pdb_names
        df_subset = filter(row -> row.File == pdb &&
            (row.Folder == pdb*"_AlphaConformer" || row.Folder == pdb * "_AF_CLUSTER" || row.Folder == pdb * "_BIOEMU_MSA"), df)

        # Prep the data for the scatter plot
        for row in eachrow(df_subset)
            if row.Folder == pdb*"_AlphaConformer"
                push!(method_labels, "AlphaConformer")
                push!(value_types, "min_apo")
                push!(values, row.min_apo)
                push!(method_labels, "AlphaConformer")
                push!(value_types, "min_holo")
                push!(values, row.min_holo)
            elseif row.Folder == pdb * "_BIOEMU_MSA"
                push!(method_labels, "BioEmu")
                push!(value_types, "min_apo")
                push!(values, row.min_apo)
                push!(method_labels, "BioEmu")
                push!(value_types, "min_holo")
                push!(values, row.min_holo)
            elseif row.Folder == pdb * "_AF_CLUSTER"
                push!(method_labels, "AF-Cluster")
                push!(value_types, "min_apo")
                push!(values, row.min_apo)
                push!(method_labels, "AF-Cluster")
                push!(value_types, "min_holo")
                push!(values, row.min_holo)
            end
        end
    end
    # gather the data
    group = [method_labels[i] * " " * value_types[i] for i in 1:length(values)]

    # create the color
    point_colors = [i % 2 == 1 ? :orange : :lightgreen for i in 1:length(values)]
    xticks_pos = [2i-1.0 for i in 1:length(group)]
    xticks_labels = ["AlphaConformer","BioEmu","AF-Cluster" ]

    p = scatter(
        group, values;
        color=point_colors,
        markerstrokewidth=0.5,
        markersize=8,
        legend=false,
        ylabel="RMSD",
        framestyle=:box,
        titlefont=font(14, "Arial"),
        tickfont=font(10),
        xticks = (xticks_pos, xticks_labels)
    )

    moyenne=[]
    for g in unique(group)
        @show g
        idx = findall(x -> x == g, group)
        @show idx
        y_mean = mean(values[idx])
        @show y_mean
        push!(moyenne, y_mean)
    end
    scatter!(group, moyenne; marker=:hline, color=:black, markersize=12, label=false)
    savefig(p, "comparaison_methodes_2_$pdb.png")
end 

"""
One plot : in abscisse the number of cluster in ordonné the performance 
One line for min apo other for min holo 
Can analyse the performance compare to the number of cluster 
"""
#Compare the number of cluster and performance 
function compare_performance_cluster(df)
    # On suppose que df.Nombre_cluster, df.min_apo, df.min_holo existent
    clusters = filter(!=(0), sort(unique(df.Nombre_cluster)))
    mean_apo = Float64[]
    mean_holo = Float64[]
    for cl in clusters
        sub = filter(row -> row.Nombre_cluster == cl, df)
        push!(mean_apo, mean(sub.min_apo))
        push!(mean_holo, mean(sub.min_holo))
    end
    @show clusters
    @show mean_apo
    @show mean_holo
    p = plot(
        clusters, mean_apo;
        label="min_apo", color=:orange, marker=:circle, lw=2,
        xlabel="Nombre de clusters", ylabel="RMSD", legend=:topright,
        title="Performance en fonction du nombre de clusters"
    )
    plot!(
        p, clusters, mean_holo;
        label="min_holo", color=:lightblue, marker=:circle, lw=2
    )

    # Afficher la moyenne sur chaque point
    scatter!(p, clusters, mean_apo; color=:orange, marker=:circle, ms=8, label=false)
    scatter!(p, clusters, mean_holo; color=:lightblue, marker=:circle, ms=8, label=false)

    savefig(p, "compare_performance_cluster.png")
end

"""
One plot : in abscisse the number of line in the MSA in ordonné the performance 
One line for min apo other for min holo 
Can analyse the performance compare to the number of line in MSA 
"""
#Compare the number of line in MSA and performance 
function compare_performance_MSA(df)
    # On suppose que df.Nombre_cluster, df.min_apo, df.min_holo existent
    clusters = filter(!=(0), sort(unique(df.Mean_MSA)))
    mean_apo = Float64[]
    mean_holo = Float64[]
    for cl in clusters
        sub = filter(row -> row.Mean_MSA == cl, df)
        push!(mean_apo, mean(sub.min_apo))
        push!(mean_holo, mean(sub.min_holo))
    end
    @show clusters
    @show mean_apo
    @show mean_holo
    p = plot(
        clusters, mean_apo;
        label="min_apo", color=:orange, marker=:circle, lw=2,
        xlabel="Mean MSA", ylabel="RMSD", legend=:topright,
        title="Performance en fonction du nombre de sequence dans les MSA"
    )
    plot!(
        p, clusters, mean_holo;
        label="min_holo", color=:lightblue, marker=:circle, lw=2
    )

    # Afficher la moyenne sur chaque point
    scatter!(p, clusters, mean_apo; color=:orange, marker=:circle, ms=8, label=false)
    scatter!(p, clusters, mean_holo; color=:lightblue, marker=:circle, ms=8, label=false)

    savefig(p, "compare_performance_MSA.png")
end

"""
One plot : in abscisse the number of sequence use in ordonné the performance 
One line for min apo other for min holo 
Can analyse the performance compare to the number of sequence
"""
#Compare the number of sequence and performance 
function compare_performance_sequence(df)
    # On suppose que df.Nombre_cluster, df.min_apo, df.min_holo existent
    clusters = filter(!=(0), sort(unique(df.Sequence_in_msa)))
    mean_apo = Float64[]
    mean_holo = Float64[]
    for cl in clusters
        sub = filter(row -> row.Sequence_in_msa == cl, df)
        push!(mean_apo, mean(sub.min_apo))
        push!(mean_holo, mean(sub.min_holo))
    end
    @show clusters
    @show mean_apo
    @show mean_holo
    p = plot(
        clusters, mean_apo;
        label="min_apo", color=:orange, marker=:circle, lw=2,
        xlabel="Sequence", ylabel="RMSD", legend=:topright,
        title="Performance en fonction du nombre de sequence"
    )
    plot!(
        p, clusters, mean_holo;
        label="min_holo", color=:lightblue, marker=:circle, lw=2
    )

    # Afficher la moyenne sur chaque point
    scatter!(p, clusters, mean_apo; color=:orange, marker=:circle, ms=8, label=false)
    scatter!(p, clusters, mean_holo; color=:lightblue, marker=:circle, ms=8, label=false)

    savefig(p, "compare_performance_sequence.png")
end

############################################################# MAIN ##################################################################
"""
Code to compare the performance of AlphaConformers 

Can use the code to compare the performance of AlphaConformers vs other pipeline like AF_cluster and BioEmu
Or can use it to compare the performance between multiple test of AlphaConformer
Use it for analyse the parameter performance 

Input :
- df : DataFrame that have the information we want to compare 
- analyse_method : Booleean to use the function for comparaison of methode or comparaison of parameter
--> If want analyse pipeline get in input the DF output by compare_output_AF.jl and analyse_method = true
--> If want analyse the parameter get in input the DF output by automatic_AlphaConformer_Test_results.jl and analyse_method=false
Output : 
Save plot in the current directory 

Run the code directly in the terminal no need to use a specific nodes 
"""

########################## Information to fill #################################
# Charger les données
df = CSV.read("comparaison_resultat_AlphaConformer.csv", DataFrame)
# Choose the comparaison you want
analyse_method = false
############################################################################

if analyse_method
    #function to compare the output of AlphaFold with different methods
    compare_methods(df)
    compare_methods_2(df)
else 
    #function to compare the output of AlphaConformers with different parameters
    compare_parameters(df)

    compare_performance_cluster(df)

    compare_performance_MSA(df)

    compare_performance_sequence(df)
end

######################################################## END #######################################################################
