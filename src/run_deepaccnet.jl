# DeepAccNet Apptainer runner
# ----------------------------

# Alphabetically-sorted cluster_* directory names directly under base_dir.
_cluster_dirs(base_dir::AbstractString) = sort(
    filter(
        name -> startswith(name, "cluster_") && isdir(joinpath(base_dir, name)),
        readdir(base_dir),
    ),
)

# Inner folder names found under one cluster directory, for the "no *.pdb found" error message.
_available_inners(cluster_dir::AbstractString) =
    sort(filter(n -> isdir(joinpath(cluster_dir, n)), readdir(cluster_dir)))

# Uses `_discover_conformers` from `src/conformer_clustering.jl` (included in
# `src/AlphaConformers.jl`) to find prediction `models/` files.

"""
    _deepaccnet_links(data_root, system) -> Vector{Tuple{String,String}}

Auto-detect the inner prediction-method folder and collect `(link_name, target_pdb)` pairs for
every `.pdb` model found under `data_root/system/cluster_*/<inner>/predictions/sequences/models/`.

# Behavior

Globs `cluster_*/*/predictions/sequences/models/*.pdb` directly under the system directory,
without assuming a fixed method name. Clusters are visited in alphabetical order; a cluster
missing the `models/` folder is skipped. Within a cluster, `*.pdb` files are visited in
alphabetical order; link names are `cluster_<N>_<model>.pdb`, so each link name is unique
regardless of visiting order. Raises a clear `ArgumentError` when `.pdb` files are found under
more than one distinct inner method folder across the system, since DeepAccNet expects a single
consistent prediction source.
"""
function _deepaccnet_links(data_root::AbstractString, system::AbstractString)
    found = try
        _discover_conformers(
            data_root,
            system;
            subdir = joinpath("cluster_*", "*", "predictions", "sequences", "models"),
        )
    catch e
        throw(ArgumentError(string(e)))
    end

    pairs = Tuple{String,String}[]
    inners = Set{String}()
    for (rel, target) in found
        parts = splitpath(rel)
        fname = parts[end]
        endswith(lowercase(fname), ".pdb") || continue
        push!(inners, parts[2])
        push!(pairs, ("$(parts[1])_$(fname)", target))
    end

    if length(inners) > 1
        throw(
            ArgumentError(
                "found .pdb predictions under more than one method folder for system '$system': " *
                join(sort(collect(inners)), ", ") *
                "; DeepAccNet needs a single consistent prediction source.",
            ),
        )
    end
    return pairs
end

"""
    _deepaccnet_command(sif_path, base_dir, links_dir, output_csv;
        n_threads=8, container_runtime="apptainer") -> Cmd

Build the container runtime invocation for the DeepAccNet `.sif` image, without executing it.

# Behavior

Binds `base_dir` (so the symlinks in `links_dir`, whose targets live under `base_dir`, resolve
inside the container), `links_dir`, and `dirname(output_csv)`; passes `links_dir` and
`output_csv` positionally per the image interface, followed by `--process N`.
"""
function _deepaccnet_command(
    sif_path::AbstractString,
    base_dir::AbstractString,
    links_dir::AbstractString,
    output_csv::AbstractString;
    n_threads::Integer = 8,
    container_runtime::String = "apptainer",
)
    runtime = _container_runtime(container_runtime)
    return `$runtime run --nv --no-home --cleanenv --bind $base_dir --bind $links_dir --bind $(dirname(output_csv)) $sif_path $links_dir $output_csv --process $n_threads`
end

"""
    _symlink_deepaccnet_pairs(links_dir, pairs)

Create `link_name => target` symlinks from `pairs` under `links_dir`, removing any stale path at
each link's location first.
"""
function _symlink_deepaccnet_pairs(
    links_dir::AbstractString,
    pairs::Vector{Tuple{String,String}},
)
    for (link_name, target) in pairs
        link_path = joinpath(links_dir, link_name)
        ispath(link_path) && rm(link_path; force = true)
        symlink(target, link_path)
    end
    return nothing
end

"""
    _run_deepaccnet(links_dir, base_dir, sif_path, output_csv, pairs;
        container_runtime="apptainer", n_threads=8, test_samples=0, dry_run=false)

Symlink `pairs` into `links_dir` and run the DeepAccNet container, writing `output_csv`.

# Behavior

1. With `test_samples > 0`, keep only the first `test_samples` pairs.
2. With `dry_run`, print the link plan and the container command and return, creating or running
   nothing.
3. Otherwise, refresh the symlinks in `links_dir` and run the container command via
   `run_cmd(...; check = true)`, then confirm `output_csv` was written.
"""
function _run_deepaccnet(
    links_dir::AbstractString,
    base_dir::AbstractString,
    sif_path::AbstractString,
    output_csv::AbstractString,
    pairs::Vector{Tuple{String,String}};
    container_runtime::String = "apptainer",
    n_threads::Integer = 8,
    test_samples::Integer = 0,
    dry_run::Bool = false,
)
    if test_samples > 0 && test_samples < length(pairs)
        pairs = pairs[1:test_samples]
    end

    cmd = _deepaccnet_command(
        sif_path,
        base_dir,
        links_dir,
        output_csv;
        n_threads = n_threads,
        container_runtime = container_runtime,
    )

    if dry_run
        println("[dry-run] would create $(length(pairs)) symlinks in $links_dir, e.g.:")
        for (link_name, target) in pairs[1:min(end, 3)]
            println("  $link_name -> $target")
        end
        println("[dry-run] would run:")
        println("  ", cmd)
        return nothing
    end

    _symlink_deepaccnet_pairs(links_dir, pairs)

    run_cmd(cmd; check = true)
    isfile(output_csv) ||
        throw(ErrorException("DeepAccNet did not write expected output CSV: $output_csv"))
    return nothing
end

"""
    run_deepaccnet(system, sif_path;
        data_root="/data/alphaconformers/pdb", out_dir="results",
        container_runtime="apptainer", n_threads=Threads.nthreads()) -> String

Run the DeepAccNet container over one AlphaConformers system and write its score table.

It flattens `data_root/system/cluster_*/<inner>/predictions/sequences/models/*.pdb` (auto-detecting
`inner`) into a flat folder of symlinks named `cluster_<N>_<model>.pdb`, runs the `.sif` image on
a GPU (`run --nv`), and returns the path to `deepaccnet_results.csv` - the score table that
`cluster_conformers(system; score_table=..., out_dir=out_dir)` consumes.

# Arguments

- `system::AbstractString`: AlphaConformers system name, matching `cluster_conformers`'s `system`.
- `sif_path::AbstractString`: Path to the DeepAccNet `.sif` Apptainer/Singularity image.

# Keywords

- `data_root::AbstractString`: Root directory holding `data_root/system/cluster_*` prediction
  folders. Defaults to `"/data/alphaconformers/pdb"`.
- `out_dir::AbstractString`: Shared AlphaConformers run directory, matching
  `cluster_conformers`'s `out_dir`. Results are written under `out_dir/system/deepaccnet/`.
  Defaults to `"results"`.
- `container_runtime::String`: Either `"apptainer"` or `"singularity"`. Defaults to
  `"apptainer"`.
- `n_threads::Int`: Number of CPUs for featurization, forwarded as `--process N`. Defaults to
  `Threads.nthreads()`.

# Returns

`output_csv = joinpath(out_dir, system, "deepaccnet", "deepaccnet_results.csv")`, the path to the
written score table.

# Throws

- `ArgumentError`: if the system directory does not exist, if it has no `cluster_*`
  subdirectories, if no `*.pdb` file is found (the available inner folders under the first
  cluster are listed in the message), or if `.pdb` files are found under more than one inner
  method folder.

# Behavior

1. Validate the system directory, the `cluster_*` subdirectories, and the auto-detected symlink
   pairs, raising a clear `ArgumentError` before any filesystem mutation.
2. Create `out_dir/system/deepaccnet/`, symlink the prediction models into a temporary directory,
   and run the container command.
3. Return `output_csv`.
"""
function run_deepaccnet(
    system::AbstractString,
    sif_path::AbstractString;
    data_root::AbstractString = "/data/alphaconformers/pdb",
    out_dir::AbstractString = "results",
    container_runtime::String = "apptainer",
    n_threads::Int = Threads.nthreads(),
)
    base_dir = joinpath(data_root, system)
    isdir(base_dir) || throw(ArgumentError("system directory not found: $base_dir"))

    clusters = _cluster_dirs(base_dir)
    isempty(clusters) && throw(ArgumentError("no cluster_* directories under $base_dir"))

    pairs = _deepaccnet_links(data_root, system)
    if isempty(pairs)
        msg = "no *.pdb prediction models found under $base_dir/cluster_*/<method>/predictions/sequences/models/"
        inners = _available_inners(joinpath(base_dir, clusters[1]))
        if !isempty(inners)
            msg *= ". Available method folders under $(clusters[1]): " * join(inners, ", ")
        end
        throw(ArgumentError(msg))
    end

    output_dir = joinpath(out_dir, system, "deepaccnet")
    output_csv = joinpath(output_dir, "deepaccnet_results.csv")
    mkpath(output_dir)

    mktempdir() do links_dir
        _run_deepaccnet(
            links_dir,
            base_dir,
            sif_path,
            output_csv,
            pairs;
            container_runtime = container_runtime,
            n_threads = n_threads,
        )
    end

    return output_csv
end
