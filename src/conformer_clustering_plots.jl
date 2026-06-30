# Conformer clustering plotting ------------------------------------------------------------
#
# Figure rendering for the conformer-clustering pipeline, kept separate from the pure algorithmic
# functions in `conformer_clustering.jl`. Every `_plot_*` helper takes precomputed data and writes
# one figure; `_render_figures` orchestrates them for one run.

# Figure styling
# -----

# Every figure is saved at this device resolution and figure size (600 dpi, 13×8 inches; Plots
# sizes in pixels at its 100 px/inch convention). The margin keeps the axis titles from being
# clipped at the figure edge (Plots has no automatic tight layout).
const _FIG_DPI = 600
const _FIG_SIZE = (1300, 800)
const _FIG_MARGIN = 10plt.mm

# Light gray used for a `-1` noise label.
const _NOISE_COLOR = plt.RGBA(0.7, 0.7, 0.7, 1.0)

# Pad a `d×N` projection up to `2×N` (filling missing rows with zeros) so the 2-D scatter plots
# always have an x and a y axis, even when the dimensionality reduction returns fewer than two
# real dimensions. Pure and deterministic.
function _ensure_2d(projection::AbstractMatrix)
    d, n = size(projection)
    d >= 2 && return projection
    out = zeros(eltype(projection), 2, n)
    d >= 1 && (out[1, :] = projection[1, :])
    return out
end

# Map cluster labels to a viridis gradient color keyed by sorted label position.
#
# `unique_labels` is any iterable of integer labels. The labels are sorted; the lowest maps to
# one end of the viridis gradient and the highest to the other (position `i/(n-1)`). A single
# label maps to the gradient midpoint. A `-1` noise label, if present, stays light gray. Returns a
# `Dict` from label to color. Pure and deterministic given the labels.
function _cluster_color_map(unique_labels)
    uniq = sort(unique(Int.(collect(unique_labels))))
    grad = [c for c in uniq if c != -1]
    n = length(grad)
    cmap = plt.cgrad(:viridis)
    color_map = Dict{Int,typeof(_NOISE_COLOR)}()
    (-1 in uniq) && (color_map[-1] = _NOISE_COLOR)
    for (i, lbl) in enumerate(grad)
        t = n <= 1 ? 0.5 : (i - 1) / (n - 1)
        color_map[lbl] = cmap[t]
    end
    return color_map
end

# A categorical palette of `k` evenly-spaced, equally-bright hues, approximating seaborn's
# "husl" palette used for the mini-cluster nodes. Returns a vector of `k` RGB colors.
function _husl_palette(k::Integer)
    return [convert(plt.RGB, plt.HSL(360 * (i - 1) / max(k, 1), 0.9, 0.65)) for i = 1:k]
end

# Graph helpers
# -----

# Minimum-spanning-tree edges of a symmetric `K×K` distance matrix (Prim's algorithm).
#
# Returns a vector of `(i, j, weight)` tuples (`i < j`), the `K-1` lightest edges that keep every
# node connected - the structural backbone used by the mini-cluster graph. Pure stdlib; no
# dependency beyond the distance matrix. An empty or singleton matrix yields no edges.
function _mst_edges(dist::AbstractMatrix)
    k = size(dist, 1)
    k < 2 && return Tuple{Int,Int,Float64}[]
    in_tree = falses(k)
    in_tree[1] = true
    edges = Tuple{Int,Int,Float64}[]
    for _ = 1:(k-1)
        best = (0, 0, Inf)
        for i = 1:k
            in_tree[i] || continue
            for j = 1:k
                in_tree[j] && continue
                if dist[i, j] < best[3]
                    best = (i, j, Float64(dist[i, j]))
                end
            end
        end
        best[1] == 0 && break
        i, j, w = best
        in_tree[j] = true
        push!(edges, (min(i, j), max(i, j), w))
    end
    return edges
end

# Average-linkage RMSD between sub-clusters of one KMeans cluster.
#
# `dist` is the member `N×N` pairwise RMSD matrix and `sub_labels` the per-member sub-cluster
# label. Returns `(uniq, mean)` where `uniq` are the sorted unique sub-labels and `mean` the
# `K×K` matrix of average between-sub-cluster RMSD (the diagonal is the mean within-cluster RMSD).
function _intercluster_rmsd(dist::AbstractMatrix, sub_labels::AbstractVector)
    uniq = sort(unique(sub_labels))
    k = length(uniq)
    mean = zeros(Float64, k, k)
    for a = 1:k, b = 1:k
        ia = findall(==(uniq[a]), sub_labels)
        ib = findall(==(uniq[b]), sub_labels)
        mean[a, b] = Statistics.mean(dist[ia, ib])
    end
    return uniq, mean
end

# Plot helpers
# -----

# Cluster-overview scatter: the PCA projection colored by KMeans cluster.
#
# Points are colored by their cluster label as a viridis gradient with a right-hand colorbar, not a
# categorical legend.
function _plot_cluster_overview(projection::AbstractMatrix, labels::AbstractVector, path)
    projection = _ensure_2d(projection)
    p = plt.scatter(
        projection[1, :],
        projection[2, :];
        marker_z = labels,
        color = :viridis,
        colorbar = true,
        colorbar_title = "KMeans cluster",
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
        xlabel = "PC 1",
        ylabel = "PC 2",
        title = "Cluster overview (PCA)",
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.savefig(p, path)
    return path
end

# Adaptive reference scatter: ref2-vs-ref1 (2 refs), embedding-vs-reference (1 ref) or the plain
# MDS embedding (0 refs), colored by KMeans cluster as a viridis gradient with a right-hand
# colorbar. `embedding` is the MDS layout used in the 0/1-reference cases; `rmsd_ref1`/`rmsd_ref2`
# carry the per-conformer reference RMSDs, and `ref_labels` names them on the axes.
function _plot_reference_scatter(
    ref_mode::Integer,
    labels::AbstractVector,
    embedding,
    rmsd_ref1,
    rmsd_ref2,
    ref_labels,
    path,
)
    common = (
        marker_z = labels,
        color = :viridis,
        colorbar = true,
        colorbar_title = "KMeans cluster",
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    if ref_mode == 2
        label1, label2 = ref_labels[1], ref_labels[2]
        p = plt.scatter(
            rmsd_ref2,
            rmsd_ref1;
            xlabel = "RMSD to $label2 reference (Å)",
            ylabel = "RMSD to $label1 reference (Å)",
            title = "$label1 vs $label2 RMSD",
            common...,
        )
    elseif ref_mode == 1
        label1 = ref_labels[1]
        embedding = _ensure_2d(embedding)
        p = plt.scatter(
            embedding[1, :],
            rmsd_ref1;
            xlabel = "Embedding dimension 1",
            ylabel = "RMSD to $label1 reference (Å)",
            title = "Embedding vs reference RMSD",
            common...,
        )
    else
        embedding = _ensure_2d(embedding)
        p = plt.scatter(
            embedding[1, :],
            embedding[2, :];
            xlabel = "MDS dimension 1",
            ylabel = "MDS dimension 2",
            title = "Conformer embedding",
            common...,
        )
    end
    plt.savefig(p, path)
    return path
end

# Score-colored RMSD scatter: the PCA projection of the scored conformers colored by score.
function _plot_rmsd_by_score(projection::AbstractMatrix, scores::AbstractVector, path)
    projection = _ensure_2d(projection)
    p = plt.scatter(
        projection[1, :],
        projection[2, :];
        marker_z = scores,
        markersize = 5,
        markerstrokewidth = 0,
        colorbar_title = "score",
        xlabel = "PC 1",
        ylabel = "PC 2",
        title = "Embedding colored by score",
        legend = false,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.savefig(p, path)
    return path
end

# KMeans clustering with DeepAccNet filtering: the cluster-overview scatter (viridis gradient)
# showing the filter result - conformers in surviving clusters keep their cluster color, those in
# filtered-out clusters are de-emphasized in light gray. `surviving` is the set of surviving
# cluster labels. Written under `deepaccnet/`.
function _plot_kmeans_deepaccnet_filtered(
    projection::AbstractMatrix,
    labels::AbstractVector,
    surviving,
    path,
)
    projection = _ensure_2d(projection)
    keep = Set(surviving)
    color_map = _cluster_color_map(unique(labels))
    point_colors = [l in keep ? color_map[l] : _NOISE_COLOR for l in labels]
    n_surv = count(in(keep), labels)
    p = plt.scatter(
        projection[1, :],
        projection[2, :];
        color = point_colors,
        markersize = 5,
        markerstrokewidth = 0,
        markeralpha = 0.7,
        label = "",
        xlabel = "PC 1",
        ylabel = "PC 2",
        title = "KMeans clusters after DeepAccNet filter ($n_surv surviving)",
        legend = false,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.savefig(p, path)
    return path
end

# Dendrogram of one cluster's agglomerative tree with a dashed line at the `threshold` cut.
#
# `tree` is a `Clustering.Hclust`, drawn with the official StatsPlots `Hclust` plot recipe. A
# `nothing` tree (a cluster too small to merge) draws a single leaf marker so the figure still
# appears.
function _plot_dendrogram(tree, threshold::Real, title::AbstractString, path)
    if tree === nothing
        p = plt.scatter(
            [1.0],
            [0.0];
            legend = false,
            xlabel = "conformer (leaf order)",
            ylabel = "merge height (Å)",
            title = title,
            size = _FIG_SIZE,
            dpi = _FIG_DPI,
            left_margin = _FIG_MARGIN,
            right_margin = _FIG_MARGIN,
            top_margin = _FIG_MARGIN,
            bottom_margin = _FIG_MARGIN,
        )
        plt.savefig(p, path)
        return path
    end

    p = plt.plot(
        tree;   
        color = :steelblue,
        legend = false,
        xlabel = "conformer (leaf order)",
        ylabel = "merge height (Å)",
        title = title,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.hline!(p, [threshold]; color = :red, linestyle = :dash)
    plt.savefig(p, path)
    return path
end

# Adaptive sub-cluster scatter: the MDS embedding of one cluster's conformers colored by their
# 1-based sub-cluster label as a viridis gradient with a right-hand colorbar.
function _plot_subcluster_scatter(
    embedding::AbstractMatrix,
    sub_labels::AbstractVector,
    title::AbstractString,
    path,
)
    embedding = _ensure_2d(embedding)
    p = plt.scatter(
        embedding[1, :],
        embedding[2, :];
        marker_z = sub_labels,
        color = :viridis,
        colorbar = true,
        colorbar_title = "sub-cluster",
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
        xlabel = "MDS dimension 1",
        ylabel = "MDS dimension 2",
        title = title,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    plt.savefig(p, path)
    return path
end

# Mini-cluster graph for one KMeans cluster: each node is an agglomerative sub-cluster, laid out
# by 2-D classical MDS of the K×K average-linkage RMSD matrix so screen distance approximates
# inter-cluster RMSD; edges are the minimum-spanning-tree backbone, colored and weighted by the
# between-sub-cluster RMSD (closer = darker/thicker) with a matching colorbar. Skipped (returns
# `nothing`) when the cluster has fewer than two sub-clusters. The node layout uses classical MDS
# (PCoA).
function _plot_minicluster_graph(
    dist::AbstractMatrix,
    sub_labels::AbstractVector,
    title::AbstractString,
    path,
)
    uniq, inter = _intercluster_rmsd(dist, sub_labels)
    k = length(uniq)
    k < 2 && return nothing

    symmetric = (inter .+ inter') ./ 2
    for i = 1:k
        symmetric[i, i] = 0.0
    end
    symmetric[symmetric .< 0] .= 0.0
    coords = _ensure_2d(_classical_mds(symmetric, 2))
    edges = _mst_edges(symmetric)
    sizes = [count(==(u), sub_labels) for u in uniq]
    palette = _husl_palette(k)        # categorical node colors, matching seaborn "husl"
    edge_cmap = plt.cgrad(:viridis; rev = true)   # viridis_r edge colormap

    weights = isempty(edges) ? [1.0] : [e[3] for e in edges]
    wmin, wmax = minimum(weights), max(maximum(weights), minimum(weights) + 1e-6)

    p = plt.plot(;
        legend = :outerbottom,
        legendtitle = "mini-cluster (N = #structures)",
        legend_column = min(k, 8, max(3, (k + 4) ÷ 5)),
        xlabel = "MDS-1 (Å)",
        ylabel = "MDS-2 (Å)",
        title = "$title\n$k mini-clusters; edge color = avg-linkage RMSD (see bar); " *
                "MST backbone ($(length(edges)) edges)",
        aspect_ratio = :equal,
        grid = true,
        gridstyle = :dash,
        gridalpha = 0.3,
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )

    # Edges: viridis_r colored and weighted by between-sub-cluster RMSD (closer = darker/thicker).
    for (i, j, w) in edges
        t = (w - wmin) / (wmax - wmin)
        plt.plot!(
            p,
            [coords[1, i], coords[1, j]],
            [coords[2, i], coords[2, j]];
            color = edge_cmap[t],
            linewidth = 0.6 + 2.2 * (1 - t),
            alpha = 0.85,
            label = "",
        )
    end
    # A hidden 2-point series carries the edge-RMSD range so a colorbar is drawn.
    if !isempty(edges)
        plt.scatter!(
            p,
            [coords[1, 1], coords[1, 1]],
            [coords[2, 1], coords[2, 1]];
            marker_z = [wmin, wmax],
            color = edge_cmap,
            markeralpha = 0,
            colorbar = true,
            colorbar_title = "edge = avg-linkage RMSD between mini-clusters (Å)",
            label = "",
        )
    end
    # Nodes: one series per mini-cluster so the legend lists each with its member count.
    for idx = 1:k
        plt.scatter!(
            p,
            [coords[1, idx]],
            [coords[2, idx]];
            color = palette[idx],
            markersize = 10,
            markerstrokewidth = 0.8,
            markerstrokecolor = :black,
            label = "$(uniq[idx]) (N=$(sizes[idx]))",
            series_annotations = plt.text(string(uniq[idx]), 8, :black),
        )
    end
    plt.savefig(p, path)
    return path
end

# Query-path view (cluster-colored): rank every KMeans cluster along x by centroid-to-centroid
# CA-RMSD to the query cluster (the query cluster is rank 0 at y = 0), with y the RMSD to the
# query cluster. Nodes are colored by the KMeans cluster gradient; reference-assigned clusters
# are overlaid with a star.
function _plot_query_path(
    inter_rmsd::AbstractMatrix,
    query_cluster::Integer,
    cluster_ids::AbstractVector,
    ref_clusters,
    path,
)
    n = length(cluster_ids)
    qpos = findfirst(==(query_cluster), cluster_ids)
    others = sort(
        [c for c in cluster_ids if c != query_cluster];
        by = c -> inter_rmsd[qpos, findfirst(==(c), cluster_ids)],
    )
    visited = vcat([query_cluster], others)
    xs = collect(0:(n-1))
    ys = [inter_rmsd[qpos, findfirst(==(c), cluster_ids)] for c in visited]

    colors = _cluster_color_map(cluster_ids)
    p = plt.plot(
        xs,
        ys;
        color = :gray,
        linewidth = 1.2,
        alpha = 0.7,
        legend = false,
        xlabel = "cluster rank by RMSD to query (0 = query cluster)",
        ylabel = "CA-RMSD to query cluster (Å)",
        title = "Query distance profile (cluster colors)",
        size = _FIG_SIZE,
        dpi = _FIG_DPI,
        left_margin = _FIG_MARGIN,
        right_margin = _FIG_MARGIN,
        top_margin = _FIG_MARGIN,
        bottom_margin = _FIG_MARGIN,
    )
    for (i, c) in enumerate(visited)
        plt.scatter!(
            p,
            [xs[i]],
            [ys[i]];
            color = colors[c],
            markersize = c == query_cluster ? 11 : 8,
            markerstrokewidth = c == query_cluster ? 2.5 : 0.8,
            markerstrokecolor = :black,
            series_annotations = plt.text(string(c), 7, :black),
        )
    end
    if ref_clusters !== nothing
        for (name, cid) in ref_clusters
            (cid === nothing || !(cid in cluster_ids)) && continue
            i = findfirst(==(cid), visited)
            plt.scatter!(
                p,
                [xs[i]],
                [ys[i]];
                color = colors[cid],
                markershape = :star5,
                markersize = 16,
                markerstrokewidth = 1.8,
                markerstrokecolor = :black,
            )
        end
    end
    plt.savefig(p, path)
    return path
end

# Figure orchestration
# -----

# Render the kept figure set for one run into the per-cluster output tree, adapting to the
# reference count and the optional score filter. Returns the vector of written figure paths.
#
# Pure orchestration over the `_plot_*` helpers; computes the PCA/MDS layouts on the fly. Sets
# `GKSwstype` so the GR backend renders offscreen on headless machines. System-level figures land
# at `system_dir`, score-filter figures under `deepaccnet_dir`, and each KMeans cluster's figures
# in its own `kmeans<C>/` folder (carried in `cluster_results`).
function _render_figures(;
    system_dir,
    deepaccnet_dir,
    ref_mode::Integer,
    labels::AbstractVector,
    features::AbstractMatrix,
    rmsd_ref1,
    rmsd_ref2,
    ref_labels,
    reference_clusters,
    scored_idx::AbstractVector,
    scores::AbstractVector,
    surviving,
    cluster_results,
    linkage::Symbol,
    threshold::Real,
    query_cluster,
    inter_rmsd,
    cluster_ids,
)
    get!(ENV, "GKSwstype", "100")
    figures = String[]

    # System-level clustering figures.
    projection = _pca(features, 2)
    push!(
        figures,
        _plot_cluster_overview(
            projection,
            labels,
            joinpath(system_dir, "cluster_overview.png"),
        ),
    )

    embedding =
        ref_mode == 2 ? zeros(2, length(labels)) :
        _classical_mds(_feature_distance_matrix(features), 2)
    push!(
        figures,
        _plot_reference_scatter(
            ref_mode,
            labels,
            embedding,
            rmsd_ref1,
            rmsd_ref2,
            ref_labels,
            joinpath(system_dir, "reference_scatter.png"),
        ),
    )

    # Query-path view - only when a query cluster was resolved (C4).
    if query_cluster !== nothing
        ref_map =
            reference_clusters === nothing ? nothing :
            Dict(
                string(r) => c for
                (r, c) in zip(reference_clusters.Reference, reference_clusters.Cluster)
            )
        push!(
            figures,
            _plot_query_path(
                inter_rmsd,
                query_cluster,
                cluster_ids,
                ref_map,
                joinpath(system_dir, "query_path.png"),
            ),
        )
    end

    # Score-filter figures (only when the filter ran).
    if deepaccnet_dir !== nothing
        push!(
            figures,
            _plot_rmsd_by_score(
                _ensure_2d(_pca(features, 2))[:, scored_idx],
                scores,
                joinpath(deepaccnet_dir, "rmsd_scatter_by_score.png"),
            ),
        )
        push!(
            figures,
            _plot_kmeans_deepaccnet_filtered(
                projection,
                labels,
                surviving,
                joinpath(deepaccnet_dir, "kmeans_deepaccnet_filtered.png"),
            ),
        )
    end

    # Per-cluster agglomerative figures, one set per `kmeans<C>/` folder.
    for cr in cluster_results
        dist = cr.dist
        sub_labels = cr.sub_labels
        cluster_dir = cr.cluster_dir
        title = "Cluster $(cr.lbl) @ $(threshold) Å cut"

        cluster_embedding = _classical_mds(dist, 2)
        push!(
            figures,
            _plot_subcluster_scatter(
                cluster_embedding,
                sub_labels,
                title,
                joinpath(cluster_dir, "subcluster_scatter.png"),
            ),
        )

        tree = size(dist, 1) >= 2 ? Clustering.hclust(dist; linkage = linkage) : nothing
        push!(
            figures,
            _plot_dendrogram(
                tree,
                threshold,
                title,
                joinpath(cluster_dir, "dendrogram.png"),
            ),
        )
        graph = _plot_minicluster_graph(
            dist,
            sub_labels,
            "Cluster $(cr.lbl) mini-clusters",
            joinpath(cluster_dir, "minicluster_graph.png"),
        )
        graph === nothing || push!(figures, graph)
    end

    return figures
end
