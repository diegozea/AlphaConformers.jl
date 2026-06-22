
"""
    safe_move(src_pattern, dst_dir)

Move all files matching a glob pattern into a destination folder.

# Arguments
- `src_pattern`: Glob pattern for files to move.
- `dst_dir`: Destination directory.

# Returns
`nothing`. Files that cannot be moved are reported and skipped.
"""
function safe_move(src_pattern::AbstractString, dst_dir::AbstractString)
    for src in glob(src_pattern)
        try
            dst = isdir(dst_dir) ? joinpath(dst_dir, basename(src)) : dst_dir
            mv(src, dst; force = true)
            println("Moved: $src → $dst_dir")
        catch e
            println("⚠️ Could not move $src: $e")
        end
    end
end

"""
    organize_files(output_dir)

Reorganize raw ColabFold output files into a structured directory hierarchy.

# Arguments
- `output_dir` : path to the folder containing raw ColabFold output files.

# Output structure
output_dir/
└── predictions/
    ├── config.json
    ├── log.txt
    ├── cite.bibtex
    └── sequences/
        ├── sequences.a3m
        ├── sequences.done.txt
        ├── models/       ← .pdb files
        ├── plots/        ← .png files
        └── scores/       ← .json files
"""
function organize_files(output_dir::AbstractString)
    println("📂 Reorganizing files in $output_dir ...")
    cwd = pwd()
    cd(output_dir)

    try
        # Main folder.
        predictions_dir = joinpath(output_dir, "predictions")
        mkdir(predictions_dir)

        mv("config.json", joinpath(predictions_dir, "config.json"))
        mv("log.txt", joinpath(predictions_dir, "log.txt"))
        mv("cite.bibtex", joinpath(predictions_dir, "cite.bibtex"))

        # Sequences folder.
        seq_dir = joinpath(predictions_dir, "sequences")
        mkdir(seq_dir)
        if isfile("sequences.a3m") || isfile("sequences.done.txt")
            mv("sequences.a3m", joinpath(seq_dir, "sequences.a3m"))
            mv("sequences.done.txt", joinpath(seq_dir, "sequences.done.txt"))
        end

        # Models folder.
        models_dir = joinpath(seq_dir, "models")
        mkdir(models_dir)

        for file in glob("*.pdb", output_dir)
            mv(file, joinpath(models_dir, basename(file)); force = true)
        end

        # Plots folder.
        plots_dir = joinpath(seq_dir, "plots")
        mkdir(plots_dir)
        for file in glob("*.png", output_dir)
            mv(file, joinpath(plots_dir, basename(file)); force = true)
        end

        # Scores folder.
        scores_dir = joinpath(seq_dir, "scores")
        mkdir(scores_dir)
        for file in glob("*.json", output_dir)
            mv(file, joinpath(scores_dir, basename(file)); force = true)
        end

        println("\n✅ Reorganization complete!")
    finally
        cd(cwd)
    end
end

"""
    organize_files_af3(output_dir, run_name)

Reorganize raw AlphaFold3 output files into a structured directory hierarchy.
Input : 
- `output_dir` : path to the folder containing raw AlphaFold3 output.
- `run_name`   : name of the run, used to locate AlphaFold3's output subfolder.
Output : 
output_dir/
└── predictions/
    ├── config.json              ← from <run_name>_data.json
    ├── cite.bibtex              ← from TERMS_OF_USE.md
    ├── ranking_scores.csv
    └── sequences/
        ├── <run_name>_model.cif
        ├── scores/
        │   ├── <run_name>_confidences.json
        │   ├── <run_name>_summary_confidences.json
        │   ├── seed_confidences.json        ← one per seed
        │   └── seed_summary_confidences.json
        └── models/
        └── seed*.cif                     ← one per seed
"""
function organize_files_af3(output_dir::AbstractString, run_name::String)
    println("📂 Reorganizing files in $output_dir ...")
    cwd = pwd()
    cd(output_dir)
    try
        output_path=joinpath(output_dir, lowercase(run_name))
        # Main folder.
        predictions_dir = joinpath(output_dir, "predictions")
        mkpath(predictions_dir)
        cp(
            joinpath(output_path, lowercase(run_name)*"_data.json"),
            joinpath(predictions_dir, "config.json");
            force = true,
        )
        cp(
            joinpath(output_path, "TERMS_OF_USE.md"),
            joinpath(predictions_dir, "cite.bibtex");
            force = true,
        )
        cp(
            joinpath(output_path, "ranking_scores.csv"),
            joinpath(predictions_dir, "ranking_scores.csv");
            force = true,
        )

        seq_dir = joinpath(predictions_dir, "sequences")
        mkpath(seq_dir)
        cp(
            joinpath(output_path, lowercase(run_name)*"_model.cif"),
            joinpath(seq_dir, lowercase(run_name)*"_model.cif");
            force = true,
        )

        scores_dir = joinpath(seq_dir, "scores")
        mkpath(scores_dir)
        cp(
            joinpath(output_path, lowercase(run_name)*"_confidences.json"),
            joinpath(scores_dir, lowercase(run_name)*"_confidences.json");
            force = true,
        )
        cp(
            joinpath(output_path, lowercase(run_name)*"_summary_confidences.json"),
            joinpath(scores_dir, lowercase(run_name)*"_summary_confidences.json");
            force = true,
        )

        models_dir = joinpath(seq_dir, "models")
        mkpath(models_dir)

        dirs=glob("seed*", output_path)
        for dir in dirs
            cp(
                joinpath(dir, "confidences.json"),
                joinpath(scores_dir, basename(dir)*"_confidences.json");
                force = true,
            )
            cp(
                joinpath(dir, "summary_confidences.json"),
                joinpath(scores_dir, basename(dir)*"_summary_confidences.json");
                force = true,
            )
            cp(
                joinpath(dir, "model.cif"),
                joinpath(models_dir, basename(dir)*".cif");
                force = true,
            )
        end

        println("\n✅ Reorganization complete!")
    finally
        cd(cwd)
    end
end
"""
    organize_files_boltz(output_dir)

Reorganize raw Boltz2 output files from multiple seeds into a structured directory hierarchy.

# Arguments
- `output_dir` : path to the folder containing `seed_*` subfolders (i.e. `bz/`).

# Output structure
output_dir/
└── predictions/
    ├── config.json       ← from processed/manifest.json (last seed)
    ├── hparams.yaml      ← from lightning_logs/version_0/hparams.yaml (last seed)
    └── sequences/
        ├── models/       ← seed_.cif
        ├── scores/       ← seed_.json
        └── plots/        ← seed_.npz
Each output file is prefixed with its seed folder name (e.g. `seed_12345_model.cif`)
to avoid collisions across seeds.
"""
function organize_files_boltz(output_dir::AbstractString)
    println("📂 Reorganizing files in $output_dir ...")
    cwd = pwd()
    cd(output_dir)

    try
        # Main folder.
        predictions_dir = joinpath(output_dir, "predictions")
        mkpath(predictions_dir)

        seq_dir = joinpath(predictions_dir, "sequences")
        mkpath(seq_dir)

        scores_dir = joinpath(seq_dir, "scores")
        mkpath(scores_dir)

        models_dir = joinpath(seq_dir, "models")
        mkpath(models_dir)

        plots_dir = joinpath(seq_dir, "plots")
        mkpath(plots_dir)

        seed_dir=glob("seed*", output_dir)
        for seed in seed_dir
            output_path=glob("*", seed)[1]
            @show output_path
            cp(
                joinpath(output_path, "processed/manifest.json"),
                joinpath(predictions_dir, "config.json");
                force = true,
            )
            cp(
                joinpath(output_path, "lightning_logs/version_0/hparams.yaml"),
                joinpath(predictions_dir, "hparams.yaml");
                force = true,
            )

            config=glob("*", joinpath(output_path, "predictions"))[1]
            models=glob("*cif", config)

            for model in models
                cp(
                    model,
                    joinpath(models_dir, basename(seed)*"_"*basename(model));
                    force = true,
                )
            end

            confidences=glob("*json", config)
            for confidence in confidences
                cp(
                    confidence,
                    joinpath(scores_dir, basename(seed)*"_"*basename(confidence));
                    force = true,
                )
            end

            plots=glob("*npz", config)
            for plot in plots
                cp(
                    plot,
                    joinpath(plots_dir, basename(seed)*"_"*basename(plot));
                    force = true,
                )
            end
        end

        println("\n✅ Reorganization complete!")
    finally
        cd(cwd)
    end
end

"""
    get_all_predictions(output_dir, folder_af2_result) -> (structures, ptm_scores)

Read all prediction models and pTM scores from cluster output folders.

# Arguments
- `output_dir`: AlphaConformers output folder containing `cluster_*` folders.
- `folder_af2_result`: Name of the prediction result folder inside each cluster.

# Returns
A pair of dictionaries. The first maps prediction names to parsed structures. The
second maps prediction names to pTM scores read from the matching JSON score files.
"""
function get_all_predictions(output_dir::String, folder_af2_result)
    clusters=glob("cluster*", output_dir)

    dic_pred_struct=Dict()
    dic_pred_ptm=Dict()
    for clu in clusters
        predictions=glob(
            "*.pdb",
            joinpath(clu, folder_af2_result, "predictions", "sequences", "models"),
        )

        for pred in predictions
            name=basename(clu)*"_"*basename(pred)
            structure=MIToS.PDB.read_file(pred, MIToS.PDB.PDBFile, group = "ATOM")
            dic_pred_struct[name]=structure
            parts = split(basename(pred), "_")
            json_name =
                "sequences_scores_" * first(splitext(join(parts[3:end], "_"))) * ".json"

            json_path = joinpath(
                clu,
                folder_af2_result,
                "predictions",
                "sequences",
                "scores",
                json_name,
            )

            if isfile(json_path)
                data = JSON3.parsefile(json_path)

                # Retrieve the pTM.
                ptm = data["ptm"]

                dic_pred_ptm[name]=ptm
            else
                @warn "Missing JSON for $pred"
            end
        end
    end

    return dic_pred_struct, dic_pred_ptm
end

"""
    compare_struct(dic_pred_struct, query_struct, cutoff_min, cutoff_max) -> Dict

Compare predicted structures with a query structure and keep models in an RMSD range.

# Arguments
- `dic_pred_struct`: Dictionary from prediction name to parsed structure.
- `query_struct`: Path to the query structure file.
- `cutoff_min`: Minimum accepted RMSD.
- `cutoff_max`: Maximum accepted RMSD.

# Returns
A dictionary mapping prediction names to RMSD values for predictions whose RMSD is
greater than `cutoff_min` and lower than `cutoff_max`.
"""
function compare_struct(
    dic_pred_struct::Dict,
    query_struct::String,
    cutoff_min::Float64,
    cutoff_max::Float64,
)
    cluster_close_query=Dict()
    @show "Read query structure"
    file_format = endswith(query_struct, ".cif") ? MIToS.PDB.MMCIFFile : MIToS.PDB.PDBFile
    query_structure=MIToS.PDB.read_file(query_struct, file_format, group = "ATOM")

    for (name, structure) in dic_pred_struct
        try
            aligned_a,
            aligned_b,
            matches,
            rmsd,
            coverage,
            identity=structural_alignment(query_structure, structure)

            if (rmsd < cutoff_max) && (rmsd > cutoff_min)

                cluster_close_query[name]=(rmsd)
            end
        catch e
            @warn "Error comparing $name to query structure: $e"
            continue
        end
    end
    @show "Number of predictions close to query: ", length(cluster_close_query)

    return cluster_close_query
end

"""
    found_uniprot_structure(search_results, sifts_uniprot_mapping, query_pdb_code,
        query_chain_code) -> DataFrame

Find Foldseek hits that belong to the same UniProt entry as a query PDB structure.

# Arguments
- `search_results`: Foldseek result table.
- `sifts_uniprot_mapping`: SIFTS UniProt mapping table.
- `query_pdb_code`: Query PDB code.
- `query_chain_code`: Query chain identifier. Currently kept for caller compatibility.

# Returns
A deduplicated subset of `search_results` whose targets start with PDB codes mapped to
the query UniProt entry.
"""
function found_uniprot_structure(
    search_results::DataFrames.DataFrame,
    sifts_uniprot_mapping::DataFrames.DataFrame,
    query_pdb_code::String,
    query_chain_code,
)
    @show query_pdb_code, query_chain_code

    query_uniprot = only(get_uniprot_acc(sifts_uniprot_mapping, query_pdb_code))

    @show query_uniprot

    query_structures = get_pdb_codes(sifts_uniprot_mapping, String(query_uniprot))

    @show query_structures

    # Columns directly.
    pdbs = String.(query_structures.PDB)
    chains = String.(query_structures.CHAIN)

    # Vectorized construction.


    uniprot_result =
        filter(r -> any(name -> startswith(r.target, name), pdbs), search_results)
    uniprot_result = unique(uniprot_result)
    @show uniprot_result
    return uniprot_result
end
"""
    found_files(file_name, possible_path1, possible_path2) -> String or nothing

Find an aligned structure file in one of two Foldseek result folders.

# Arguments
- `file_name`: Target name to search for.
- `possible_path1`: First result folder to inspect.
- `possible_path2`: Second result folder to inspect.

# Returns
The first matching PDB path, or `nothing` if no matching file is found.
"""
function found_files(file_name, possible_path1, possible_path2)
    files1=glob("*.pdb", joinpath(possible_path1, "aligned_structures"))
    file_found = filter(f -> occursin(file_name, f), files1)
    if !isempty(file_found)
        return file_found[1]
    end
    files2=glob("*.pdb", joinpath(possible_path2, "aligned_structures"))
    file_found = filter(f -> occursin(file_name, f), files2)
    if !isempty(file_found)
        return file_found[1]
    end
    return nothing
end
"""
    compare_alternative_structures(search_results, sifts_uniprot_mapping, output_dir)
    -> (min_rmsd, max_rmsd) or nothing

Estimate the RMSD range among known related structures in Foldseek results.

# Arguments
- `search_results`: Foldseek result table containing related structures.
- `sifts_uniprot_mapping`: SIFTS UniProt mapping table.
- `output_dir`: AlphaConformers output folder containing Foldseek result folders.

# Returns
A tuple with the minimum and maximum per-entry RMSD ranges. Returns `nothing` when no
related structures can be compared.
"""
function compare_alternative_structures(search_results, sifts_uniprot_mapping, output_dir)
    range_rmsd=[]
    file_analysed=Set{String}()
    for result in eachrow(search_results)
        file_name=result.target

        if file_name in file_analysed
            continue
        end
        push!(file_analysed, file_name)
        parts=split(file_name, "_")
        query_pdb_code=String(first(splitext(parts[1])))
        if length(parts)==2
            query_chain_code=String(parts[end])
        else
            query_chain_code=missing
        end
        center_uniprot=found_uniprot_structure(
            search_results,
            sifts_uniprot_mapping,
            query_pdb_code,
            query_chain_code,
        )
        @show "Number of alternative structures found for $file_name: ",
        nrow(center_uniprot)
        if isempty(center_uniprot)
            continue
        end
        max_rmsd=0.0
        for i = 1:nrow(center_uniprot)
            path_i=found_files(
                center_uniprot[i, :target],
                joinpath(output_dir, "target_db_results"),
                joinpath(output_dir, "fullpdb_mmcif_files_results"),
            )
            push!(file_analysed, center_uniprot[i, :target])
            struct1 = MIToS.PDB.read_file(path_i, MIToS.PDB.PDBFile, group = "ATOM")
            for j = (i+1):nrow(center_uniprot)
                path_j=found_files(
                    center_uniprot[j, :target],
                    joinpath(output_dir, "target_db_results"),
                    joinpath(output_dir, "fullpdb_mmcif_files_results"),
                )
                push!(file_analysed, center_uniprot[j, :target])
                struct2=MIToS.PDB.read_file(path_j, MIToS.PDB.PDBFile, group = "ATOM")

                aligned_a,
                aligned_b,
                matches,
                rmsd,
                coverage,
                identity=structural_alignment(struct1, struct2)
                if rmsd > max_rmsd
                    @show "RMSD between ",
                    center_uniprot[i, :target],
                    " and ",
                    center_uniprot[j, :target],
                    " : ",
                    rmsd
                    max_rmsd = rmsd
                end
            end
        end
        push!(range_rmsd, max_rmsd)
    end
    if isempty(range_rmsd)
        return nothing
    end
    @show "RMSD range between alternative structures: ", range_rmsd
    return minimum(range_rmsd), maximum(range_rmsd)
end

"""
    found_best_prediction(output_dir, query_struct, sifts_uniprot_mapping, folder_af2_result)
    -> (cluster_close_objectif, dic_pred_ptm) 

Identify the best structure predictions by comparing them to known experimental structures
retrieved via Foldseek, and filtering to those structurally close to the query.
Input : 
- `output_dir`            : working directory containing Foldseek results and prediction folders.
- `query_struct`          : path to the query structure (PDB file) used as structural reference.
- `sifts_uniprot_mapping` : SIFTS mapping table linking UniProt accessions to PDB entries.
- `folder_af2_result`     : subfolder name containing the AlphaFold2/ColabFold predictions.
Behavior
1. Loads all predictions and their pTM scores from `folder_af2_result`.
2. Reads Foldseek hits from `fullpdb_mmcif_files_results/*.m8` and optionally
   from `target_db_results/*.m8`, merging known UniProt structures from both.
3. Computes the RMSD range of known alternative structures via `compare_alternative_structures`.
4. Filters predictions to those within `(max_rmsd + 1)` Å of the query structure.
5. Returns only predictions that pass the RMSD threshold, along with their pTM scores.
Output : 
A tuple `(cluster_close_objectif, dic_pred_ptm)` where:
- `cluster_close_objectif` : dict mapping prediction names to their RMSD to the query.
- `dic_pred_ptm`           : dict mapping the same prediction names to their pTM scores.
Both dicts are empty if no predictions are found, no alternative structures are available,
or no predictions fall within the RMSD threshold.
"""
function found_best_prediction(
    output_dir::String,
    query_struct::String,
    sifts_uniprot_mapping,
    folder_af2_result,
)
    dic_pred_struct, dic_pred_ptm=get_all_predictions(output_dir, folder_af2_result)

    if isempty(dic_pred_struct)
        @warn "No predictions found in $output_dir for folder $folder_af2_result"
        return Dict(), Dict()
    end

    full_m8file=glob("*.m8", joinpath(output_dir, "fullpdb_mmcif_files_results"))
    @show "Reading foldseek all_search results PDB from ", full_m8file[1]
    all_search_results=read_foldseek_search_results(full_m8file[1])
    all_alternative_files=known_uniprot_structures(
        sifts_uniprot_mapping,
        all_search_results,
    )
    if isdir(joinpath(output_dir, "target_db_results"))
        @show "Comparing to known structures in target_db_results folder"
        m8file=glob("*.m8", joinpath(output_dir, "target_db_results"))
        if !isempty(m8file) && isfile(m8file[1])
            @show "Reading foldseek search results from ", m8file[1]
            search_results=read_foldseek_search_results(m8file[1])

            alternative_files=known_uniprot_structures(
                sifts_uniprot_mapping,
                search_results,
            )

            all_alternative_files=union(all_alternative_files, alternative_files)
        end
    end
    if isempty(all_alternative_files)
        @warn "No alternative structures found in foldseek results. Skipping foldseek comparison."
        return Dict(), Dict()
    end
    for fname in alternative_files

        filtered_results=filter(r -> r.target == fname, all_search_results)
        append!(search_results, filtered_results)
    end

    if nrow(search_results) < 2
        @warn "No alternative structures found in foldseek results. Skipping foldseek comparison."
        return Dict(), Dict()
    end
    @show "Comparing alternative structures to query structure"
    @show "Number of alternative structures found: ", nrow(search_results)
    output =
        compare_alternative_structures(search_results, sifts_uniprot_mapping, output_dir)
    if isempty(output)
        return Dict(), Dict()
    else
        min_rmsd, max_rmsd = output
    end
    cluster_close_objectif=compare_struct(dic_pred_struct, query_struct, 0.0, (max_rmsd+1))
    if isempty(cluster_close_objectif)
        @warn "No predictions close to query within $(max_rmsd+1). Returning empty results."
        return Dict(), Dict()
    end
    keyset = Set(keys(cluster_close_objectif))
    dic_pred_ptm_close_objectif_all = Dict(k => v for (k, v) in dic_pred_ptm if k in keyset)

    return cluster_close_objectif, dic_pred_ptm_close_objectif_all

end

"""
    triage_outputs(output_dir, query_struct, sifts_uniprot_mapping, folder_af2_result)
    -> (selected_predictions, ptm_scores)

Triage prediction outputs from an AlphaConformers run.

This is the official final step of the pipeline. For now, it calls
`found_best_prediction` and returns the same result. Future versions may add more
triage logic on top of that filtering step.

# Arguments
- `output_dir`: AlphaConformers output folder containing Foldseek results and
  prediction folders.
- `query_struct`: Query structure used as the structural reference.
- `sifts_uniprot_mapping`: SIFTS mapping table linking UniProt accessions to PDB
  entries.
- `folder_af2_result`: Prediction result folder inside each cluster.

# Returns
A tuple with the selected predictions and their pTM scores.
"""
function triage_outputs(
    output_dir::String,
    query_struct::String,
    sifts_uniprot_mapping,
    folder_af2_result,
)
    return found_best_prediction(
        output_dir,
        query_struct,
        sifts_uniprot_mapping,
        folder_af2_result,
    )
end
