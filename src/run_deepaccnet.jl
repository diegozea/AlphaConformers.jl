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
    _deepaccnet_links(base_dir, inner) -> Vector{Tuple{String,String}}

Collect `(link_name, target_pdb)` pairs for every prediction model found under
`base_dir/cluster_*/<inner>/predictions/sequences/models/`.

# Behavior

Clusters are visited in alphabetical order; a cluster missing the `models/` folder for `inner`
is skipped. Within a cluster, `*.pdb` files are visited in alphabetical order; link names are
`cluster_<N>_<model>.pdb`, so each link name is unique regardless of visiting order.
"""
function _deepaccnet_links(base_dir::AbstractString, inner::AbstractString)
    # base_dir is the system directory; use its parent as data_root and its basename as system
    data_root = dirname(base_dir)
    system = basename(base_dir)
    found = try
        _discover_conformers(data_root, system; method = inner)
    catch e
        throw(ArgumentError(string(e)))
    end
    
    pairs = Tuple{String,String}[]
    for (rel, abs) in found
        parts = splitpath(rel)
        cluster = parts[1]
        fname = parts[end]
        push!(pairs, ("$(cluster)_$(fname)", abs))
    end
    return pairs
end

"""
    _deepaccnet_command(sif_path, base_dir, links_dir, output_csv; process=8) -> Cmd

Build the `apptainer run --nv` invocation for the DeepAccNet `.sif` image, without executing it.

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
    process::Integer = 8,
)
    return `apptainer run --nv --no-home --cleanenv --bind $base_dir --bind $links_dir --bind $(dirname(output_csv)) $sif_path $links_dir $output_csv --process $process`
end

"""
    _clean_deepaccnet_intermediates(links_dir) -> Int

Delete the heavy DeepAccNet intermediates (`bert_*.npy`, `*.features.npz`, `*.fa`) from
`links_dir`. Returns the number of files removed.
"""
function _clean_deepaccnet_intermediates(links_dir::AbstractString)
    removed = 0
    for f in readdir(links_dir)
        if (startswith(f, "bert_") && endswith(f, ".npy")) ||
           endswith(f, ".features.npz") ||
           endswith(f, ".fa")
            rm(joinpath(links_dir, f); force = true)
            removed += 1
        end
    end
    return removed
end

"""
    _gpu_available() -> Bool

Return whether a GPU usable by `apptainer run --nv` is present, by checking that `nvidia-smi`
is on `PATH` and runs successfully. Never throws; returns `false` on any failure.
"""
function _gpu_available()
    try
        Sys.which("nvidia-smi") === nothing && return false
        return success(`nvidia-smi`)
    catch
        return false
    end
end

"""
    run_deepaccnet(base_dir, sif_path;
        inner="af", output_dir=joinpath(base_dir, "deepaccnet_output"),
        links_dir=joinpath(output_dir, "links"), process=8, test_samples=0,
        require_gpu=true, keep_intermediates=false, dry_run=false) -> String

Run the DeepAccNet Apptainer image over one AlphaConformers system and write its score table.

It flattens `base_dir/cluster_*/<inner>/predictions/sequences/models/*.pdb` into a flat folder of
symlinks named `cluster_<N>_<model>.pdb`, runs the `.sif` image on a GPU
(`apptainer run --nv`), and returns the path to `deepaccnet_results.csv` - the score table that
`cluster_conformers(system; score_table=...)` consumes.

# Arguments

- `base_dir::AbstractString`: One AlphaConformers system directory containing `cluster_*`
  subfolders.
- `sif_path::AbstractString`: Path to the DeepAccNet `.sif` Apptainer image.

# Keywords

- `inner::AbstractString`: Inner folder name under each `cluster_*` holding the prediction
  models (e.g. `"af_pdb"`, `"af_cif"`). Defaults to `"af"`.
- `output_dir::AbstractString`: Directory for the final `deepaccnet_results.csv`. Defaults to
  `joinpath(base_dir, "deepaccnet_output")`.
- `links_dir::AbstractString`: Flat symlink and intermediate directory. Defaults to
  `joinpath(output_dir, "links")`.
- `process::Integer`: Number of CPUs for featurization, forwarded as `--process N`. Defaults
  to `8`.
- `test_samples::Integer`: Limit to the first N symlink pairs, for a quick smoke test. `0`
  (the default) uses all of them.
- `require_gpu::Bool`: Whether to refuse a real run when no GPU is detected. Defaults to `true`.
- `keep_intermediates::Bool`: Whether to keep `bert_*.npy` / `*.features.npz` / `*.fa` in
  `links_dir` after a successful run. Defaults to `false`.
- `dry_run::Bool`: Print the link plan and the apptainer command without creating or running
  anything; bypasses the GPU check. Defaults to `false`.

# Returns

`output_csv = joinpath(output_dir, "deepaccnet_results.csv")`, the path to the written (or, in
`dry_run` mode, the would-be) score table.

# Throws

- `ArgumentError`: if `base_dir` does not exist, if it has no `cluster_*` subdirectories, if no
  `*.pdb` file is found for `inner` (the available inner folders under the first cluster are
  listed in the message), or - for a real run with `require_gpu=true` - if no GPU is detected.

# Behavior

1. Validate `base_dir`, the `cluster_*` subdirectories, and the symlink pairs for `inner`,
   raising a clear `ArgumentError` before any filesystem mutation.
2. With `dry_run`, print the link plan and the apptainer command and return `output_csv`,
   bypassing the GPU check and creating/running nothing.
3. Otherwise, check `require_gpu` against `_gpu_available()`, refusing with an `ArgumentError`
   before any filesystem mutation when no GPU is available.
4. Create `links_dir`, refresh the symlinks (removing any stale path first), create
   `output_dir`, and run the apptainer command via `run_cmd`.
5. Unless `keep_intermediates`, delete the heavy intermediates from `links_dir`.
6. Return `output_csv`.
"""
function run_deepaccnet(
    base_dir::AbstractString,
    sif_path::AbstractString;
    inner::AbstractString = "af",
    output_dir::AbstractString = joinpath(base_dir, "deepaccnet_output"),
    links_dir::AbstractString = joinpath(output_dir, "links"),
    process::Integer = 8,
    test_samples::Integer = 0,
    require_gpu::Bool = true,
    keep_intermediates::Bool = false,
    dry_run::Bool = false,
)
    isdir(base_dir) || throw(ArgumentError("base_dir not found: $base_dir"))

    clusters = _cluster_dirs(base_dir)
    isempty(clusters) && throw(ArgumentError("no cluster_* directories under $base_dir"))

    pairs = _deepaccnet_links(base_dir, inner)
    if isempty(pairs)
        msg = "no *.pdb found for inner='$inner' under $base_dir"
        inners = _available_inners(joinpath(base_dir, clusters[1]))
        if !isempty(inners)
            msg *= ". Available inner folders under $(clusters[1]): " * join(inners, ", ")
        end
        throw(ArgumentError(msg))
    end

    if test_samples > 0 && test_samples < length(pairs)
        pairs = pairs[1:test_samples]
    end

    output_csv = joinpath(output_dir, "deepaccnet_results.csv")

    if dry_run
        println("[dry-run] would create $(length(pairs)) symlinks in $links_dir, e.g.:")
        for (link_name, target) in pairs[1:min(end, 3)]
            println("  $link_name -> $target")
        end
        cmd = _deepaccnet_command(
            sif_path,
            base_dir,
            links_dir,
            output_csv;
            process = process,
        )
        println("[dry-run] would run:")
        println("  ", cmd)
        return output_csv
    end

    if require_gpu && !_gpu_available()
        throw(
            ArgumentError(
                "DeepAccNet needs apptainer run --nv (no GPU detected); pass require_gpu=false to override",
            ),
        )
    end

    mkpath(links_dir)
    for (link_name, target) in pairs
        link_path = joinpath(links_dir, link_name)
        ispath(link_path) && rm(link_path; force = true)
        symlink(target, link_path)
    end
    mkpath(output_dir)

    run_cmd(
        _deepaccnet_command(sif_path, base_dir, links_dir, output_csv; process = process),
    )

    keep_intermediates || _clean_deepaccnet_intermediates(links_dir)

    return output_csv
end
