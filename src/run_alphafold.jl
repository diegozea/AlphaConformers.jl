

# === FONCTIONS UTILITAIRES ===
"""
    run_cmd(cmd)

Run a shell command, logging a success message or a warning on failure.
Input : 
- `cmd` : a Julia `Cmd` object to execute.
Errors are caught and logged as warnings so the calling code can continue.
"""
function run_cmd(cmd::Cmd)
    try
        run(cmd)
        println("✅ end without error\n")
    catch e
        @warn "⚠️ Error of execution for : $e"
    end
end

"""
    run_alphafold_one_run(clusters_folder, SIF_PATH, CACHE_DIR)

Run ColabFold structure predictions for all `cluster_*` subfolders, using
custom structural templates from each cluster's `templates_adaptative/` folder.
Run with .sif image for ColabFold 1.5.5
Input : 
- `clusters_folder` : path to the folder containing `cluster_*` subfolders.
- `SIF_PATH`        : path to the ColabFold Apptainer `.sif` image.
- `CACHE_DIR`       : path to the cache directory for Apptainer.
"""
function run_alphafold_one_run(clusters_folder::String, SIF_PATH::String, CACHE_DIR::String)
    ### Check input directories
    if isempty(clusters_folder)
        throw(ErrorException("No cluster_* folders were found in $clusters_folder"))
    end
    input_dirs = sort(
        filter(
            d -> isdir(d) && startswith(basename(d), "cluster_"),
            readdir(clusters_folder; join = true),
        ),
    )

    println("📂 Number of folder found: ", length(input_dirs))
    seed = rand(10_000:99_999)
    #seed = 96826
    for input_dir in input_dirs
        #Check template info 

        name = basename(input_dir)
        output_dir = joinpath(input_dir, "af")
        @show "Output directory: $output_dir"
        if isdir(output_dir)
            if isdir(joinpath(output_dir, "predictions", "sequences", "models"))
                @show "Folder already process"
                continue
            else
                rm(output_dir; recursive = true, force = true)
            end
        end
        mkdir(output_dir)

        println("🚀Start $name ...")

        cmd = `apptainer exec --nv --no-home --cleanenv \
            --bind $input_dir:/mnt/input \
            --bind $output_dir:/mnt/output \
            --bind $CACHE_DIR:/cache \
            $SIF_PATH \
            colabfold_batch \
                /mnt/input/sequences.a3m /mnt/output \
                --templates \
                --custom-template-path /mnt/input/templates_adaptative/ \
                --random-seed $seed \
                --num-seeds 5 \
                --use-dropout \
                --num-models 2 \
                --overwrite-existing-results`
        run_cmd(cmd)

        organize_files(output_dir)

    end

    println("🎉 All the ColabFold run are finish !")
end

"""
    get_msa_sequence_afdb(uniprot_id, output_dir) -> path or nothing

Download the precomputed MSA for a UniProt entry from the AlphaFold Database.
Input : 
- `uniprot_id` : UniProt accession (e.g. `"P00533"`).
- `output_dir` : folder where the `.a3m` file will be saved.
Output : 
The path to the saved `<UNIPROT_ID>_msa.a3m` file, or `nothing` if failed
If multiple AlphaFold entries exist for the same UniProt ID, the first one is used.
"""
function get_msa_sequence_afdb(uniprot_id::String, output_dir::String)

    json_result = nothing
    try
        json_result = MIToS.PDB.query_alphafolddb(uniprot_id)
    catch e
        if occursin("multiple elements", sprint(showerror, e))
            println(
                "⚠️ Multiple AlphaFold entries found for $uniprot_id. Using the first one.",
            )
            resp = HTTP.get("https://alphafold.ebi.ac.uk/api/prediction/$uniprot_id")
            data = JSON3.read(resp.body)
            json_result = data[1]
        else
            @warn "Error querying AlphaFold DB for $uniprot_id : ", e
            return nothing
        end
    end
    if json_result === nothing
        @warn "No AlphaFold data available for $uniprot_id"
        return nothing
    end
    if !haskey(json_result, "msaUrl")
        @warn "No msaUrl for $uniprot_id"
        return nothing
    end
    msa_url = json_result["msaUrl"]

    msa_path_save=joinpath(output_dir, "$(uppercase(uniprot_id))_msa.a3m")
    msa_file = MIToS.Utils.download_file(msa_url, msa_path_save)
    return msa_path_save

end

"""
    parse_plddt_info_in_line(line) -> NamedTuple or nothing

Parse one ColabFold log line containing ranking and confidence scores.

# Arguments
- `line`: Log line to inspect.

# Returns
A named tuple with `rank`, `pLDDT`, and `pTM` when the line contains all three
fields. Returns `nothing` for unrelated lines.
"""
function parse_plddt_info_in_line(line::String)
    @show line
    m = match(r"(rank_[^ ]+).*pLDDT=(\d+\.\d+).*pTM=(\d+\.\d+)", line)
    @show m
    if m !== nothing
        return (
            rank = m.captures[1],
            pLDDT = parse(Float64, m.captures[2]),
            pTM = parse(Float64, m.captures[3]),
        )
    else
        return nothing
    end
end

"""
    run_alphafold_input_structure(uniprot, output_path, SIF_PATH, CACHE_DIR; msa_path_save)
    -> path or nothing

Run a ColabFold structure prediction for a UniProt entry using its AlphaFold DB MSA,
and return the path to the best predicted model.
Input : 
- `uniprot`     : UniProt accession (e.g. `"P00533"`).
- `output_path` : working directory where inputs and outputs are stored.
- `SIF_PATH`    : path to the ColabFold Apptainer `.sif` image.
- `CACHE_DIR`   : path to the cache directory for Apptainer.
- `msa_path_save` : if provided, use this `.a3m` file instead of downloading from AlphaFold DB.
                    The file is copied into `output_path` before the run.

Downloads (or copies) the MSA to `output_path`
Runs ColabFold with 5 models × 5 seeds and dropout enabled.
Reorganizes output files with `organize_files`.
Parses `log.txt` to extract pLDDT scores, saves them to `scores.csv`,
   and copies the best-ranked model to `output_path/best_model.pdb`.
Output : 
Path to `best_model.pdb`, or `nothing` if the folder was already processed.
"""
function run_alphafold_input_structure(
    uniprot::String,
    output_path::String,
    SIF_PATH::String,
    CACHE_DIR::String;
    msa_path_save::Union{String,Missing} = missing,
)
    ## Get afdb a3M 
    if ismissing(msa_path_save)
        msa_path_save=get_msa_sequence_afdb(uniprot, output_path)
        @show "Saved MSA path: $msa_path_save"
    else
        cp(msa_path_save, joinpath(output_path, uniprot*"_msa.a3m"))
        msa_path_save=joinpath(output_path, uniprot*"_msa.a3m")
    end
    ## Run colabfold with the msa 
    output_dir = joinpath(output_path, "af_input")
    @show "Output directory: $output_dir"
    if isdir(output_dir)
        if isdir(joinpath(output_dir, "predictions", "sequences", "models"))
            @show "Folder already process"
            return
        else
            rm(output_dir; recursive = true, force = true)
        end
    end
    mkdir(output_dir)

    seed = rand(10_000:99_999)
    cmd = `apptainer exec --nv --no-home --cleanenv \
        --bind $output_path:/mnt/input \
        --bind $output_dir:/mnt/output \
        --bind $CACHE_DIR:/cache \
        $SIF_PATH \
        colabfold_batch \
            /mnt/input/$(basename(msa_path_save)) /mnt/output \
            --random-seed $seed \
            --num-seeds 5 \
            --use-dropout \
            --num-models 5 \
            --overwrite-existing-results`

    run_cmd(cmd)

    organize_files(output_dir)

    ## Create score file 
    log_file = joinpath(output_dir, "predictions", "log.txt")

    parsed_data = []
    for line in readlines(log_file)
        if occursin("rank_", line)
            parsed = parse_plddt_info_in_line(line)
            @show parsed
            if parsed !== nothing
                push!(parsed_data, parsed)
            end
        end
    end

    df = DataFrame(parsed_data)
    @show df
    # Sort by descending pLDDT.
    sort!(df, :pLDDT, rev = true)

    CSV.write(joinpath(output_dir, "predictions", "scores.csv"), df)

    ## Select the best model
    best_model = df[1, :rank]
    file_name="$(uniprot)_msa_unrelaxed_$(best_model).pdb"
    src = joinpath(output_dir, "predictions", "sequences", "models", file_name)
    dst = joinpath(output_path, "best_model.pdb")

    cp(src, dst; force = true)

    return dst
end

### AlphaFold3 function ###

"""
    write_af3_json(out_json; run_name, sequence, chain_id="A", templates_info, seed)
    -> out_json

Write an AlphaFold3 JSON input file for one protein sequence.

# Arguments
- `out_json`: Path where the JSON file will be written.

# Keywords
- `run_name`: Name used by AlphaFold3 for this run.
- `sequence`: Protein sequence for the query chain.
- `chain_id`: Query chain identifier.
- `templates_info`: Template records in the format expected by AlphaFold3.
- `seed`: First model seed. Four following seeds are added automatically.

# Returns
The path passed as `out_json`.
"""
function write_af3_json(
    out_json::String;
    run_name::String,
    sequence::String,
    chain_id::String = "A",
    templates_info::Vector{Any},
    seed::Int = rand(10_000:99_999),
)
    L = length(sequence)

    config = Dict(
        "name" => run_name,  # MUST be unique, no extensions
        "modelSeeds" => [seed, seed+1, seed+2, seed+3, seed+4],
        "sequences" => [
            Dict(
                "protein" => Dict(
                    "id" => chain_id,
                    "sequence" => sequence,
                    "unpairedMsaPath" => "/root/af_input/sequences.a3m",
                    "pairedMsa" => "",
                    "templates" => templates_info,
                ),
            ),
        ],
        "dialect" => "alphafold3",
        "version" => 2,
    )

    open(out_json, "w") do io
        JSON3.write(io, config; indent = 2)
    end

    return out_json
end

"""
    aligned_indices(query_aln, template_aln) -> (query_indices, template_indices)

Return zero-based residue indices for columns where two aligned sequences both have
residues.

# Arguments
- `query_aln`: Aligned query sequence, using `-` for gaps.
- `template_aln`: Aligned template sequence, using `-` for gaps.

# Returns
Two integer vectors with matching query and template residue indices. The indices are
zero-based because AlphaFold3 template JSON uses zero-based indexing.
"""
function aligned_indices(query_aln::String, template_aln::String)
    @assert length(query_aln) == length(template_aln)

    query_indices = Int[]
    template_indices = Int[]

    qi = 0  # query index (0-based)
    ti = 0  # template index (0-based)

    for i in eachindex(query_aln)
        q = query_aln[i]
        t = template_aln[i]

        if q != '-' && t != '-'
            push!(query_indices, qi)
            push!(template_indices, ti)
        end

        if q != '-'
            qi += 1
        end
        if t != '-'
            ti += 1
        end
    end

    return query_indices, template_indices
end

"""
    run_alphafold3(clusters_folder, SIF_IMAGE_PATH, MODEL_PARAMETERS_DIR, DB_DIR)

Run AlphaFold3 predictions for all `cluster_*` subfolders in a given directory.
Tested with version 22b9ab8
Input : 
- `clusters_folder`    : path to the folder containing `cluster_*` subfolders.
- `SIF_IMAGE_PATH`     : path to the folder containing the `alphafold3.sif` Apptainer image.
- `MODEL_PARAMETERS_DIR` : path to the AlphaFold3 model weights.
- `DB_DIR`             : path to the AlphaFold3 genetic databases.
Output : 
For each `cluster_*` subfolder:
1. Patches any `.cif` template files for AlphaFold3 compatibility.
2. Skips folders that have already been processed (i.e. `af3/predictions/sequences/models/` exists).
3. Reads the query sequence and MSA from `sequences.a3m`.
4. Writes an `af3_config.json` input file with structural templates if a `templates/`
   subfolder is present, without otherwise.
5. Runs AlphaFold3 via Apptainer and reorganizes the output with `organize_files_af3`.
"""
function run_alphafold3(
    clusters_folder::String,
    SIF_IMAGE_PATH::String,
    MODEL_PARAMETERS_DIR::String,
    DB_DIR::String,
)
    if isempty(clusters_folder)
        throw(ErrorException("No cluster_* folders were found in $clusters_folder"))
    end

    input_dirs = sort(
        filter(
            d -> isdir(d) && startswith(basename(d), "cluster_"),
            readdir(clusters_folder; join = true),
        ),
    )
    #=
    input_dirs = sort(
        filter(
            d -> isdir(d),
            readdir(clusters_folder; join=true)
        )
    )
    =#
    println("📂 Number of folder found: ", length(input_dirs))
    run_name = split(basename(clusters_folder), "_")[1]
    seed = rand(10_000:99_999)
    for input_dir in input_dirs

        output_dir=joinpath(input_dir, "af3")
        #Get the input informations
        template_names=glob("*.cif", joinpath(input_dir, "templates"))
        for template_name in template_names
            if !startswith(basename(template_name), "t00")
                @show "Change template name for AlphaFold compatibility $(basename(template_name))"
                patch_mmcif_for_alphafold(template_name, template_name)
                rm(output_dir; recursive = true, force = true)
            end

        end

        #Create output folder
        output_dir=joinpath(input_dir, "af3")
        if isdir(output_dir)
            if isdir(joinpath(output_dir, "predictions", "sequences", "models"))
                @show "Folder already process"
                continue
            else
                rm(output_dir; recursive = true, force = true)
            end
        end
        mkdir(output_dir)

        cluster_name = basename(input_dir)
        out_json=joinpath(input_dir, "af3_config.json")

        ref_msa_path = joinpath(input_dir, "sequences.a3m")
        ids, sequences = read_a3m(ref_msa_path)
        query_id = ids[1]
        chain_id = split(query_id, "_")[2]
        sequence = sequences[1]
        templates_info = []
        if isdir(joinpath(input_dir, "templates"))
            template_names=glob("*.cif", joinpath(input_dir, "templates"))
            for template_name in template_names
                template_path =
                    joinpath("/root/af_input/templates", basename(template_name))
                split_name=String(first(splitext(basename(template_name))))
                @show split_name
                matches = [
                    sequences[i] for
                    i in eachindex(ids) if startswith(ids[i], lowercase(split_name))
                ]
                if isempty(matches)
                    error("No template found for $split_name")
                end
                template_seq = matches[1]
                real_name = [
                    ids[i] for
                    i in eachindex(ids) if startswith(ids[i], lowercase(split_name))
                ]
                chain_template=split(real_name[1], "_")
                if length(chain_template) > 1
                    chain_template = String(uppercase(chain_template[2]))
                else
                    chain_template = "A"
                end
                query_indices, template_indices = aligned_indices(sequence, template_seq)
                push!(
                    templates_info,
                    Dict(
                        "mmcifPath" => template_path,
                        "chainId" => chain_template,
                        "queryIndices" => query_indices,
                        "templateIndices" => template_indices,
                    ),
                )
            end
            # Right input file
            write_af3_json(
                out_json;
                run_name = String(run_name),
                sequence = sequence,
                chain_id = String(chain_id),
                templates_info = templates_info,
                seed = seed,
            )
        else
            # Right input file without template
            config = Dict(
                "name" => run_name,  # MUST be unique, no extensions
                "modelSeeds" => [seed, seed+1, seed+2, seed+3, seed+4],
                "sequences" => [
                    Dict(
                        "protein" => Dict(
                            "id" => chain_id,
                            "sequence" => sequence,
                            "unpairedMsaPath" => "/root/af_input/sequences.a3m",
                            "pairedMsa" => "",
                            "templates" => [],
                        ),
                    ),
                ],
                "dialect" => "alphafold3",
                "version" => 2,
            )

            open(out_json, "w") do io
                JSON3.write(io, config; indent = 2)
            end
        end



        println("🚀Start $cluster_name ...")
        # Run AlphaFold3
        cmd = `apptainer exec \
            --nv --no-home --cleanenv \
            --bind $input_dir:/root/af_input \
            --bind $output_dir:/root/af_output \
            --bind $MODEL_PARAMETERS_DIR:/root/models \
            --bind $DB_DIR:/root/public_databases \
            $SIF_IMAGE_PATH/alphafold3.sif \
            python /app/alphafold/run_alphafold.py \
            --json_path=/root/af_input/af3_config.json \
            --model_dir=/root/models \
            --db_dir=/root/public_databases \
            --output_dir=/root/af_output
                `
        try
            run(cmd)
            @info "✅ AlphaFold3 run finished for $cluster_name"
            #organize output files
            organize_files_af3(output_dir, String(run_name))
        catch e
            @warn e
        end
    end
end

#Run Boltz
"""
    run_boltz2(clusters_folder, SIF_IMAGE_PATH, CACHE_DIR_Boltz)

Run Boltz2 structure predictions for all subfolders in a given directory.
Tested with version fd69c81
Input : 
- `clusters_folder`  : path to the folder containing cluster subfolders.
- `SIF_IMAGE_PATH`   : path to the Boltz2 Apptainer `.sif` image.
- `CACHE_DIR_Boltz`  : path to the cache directory (Numba, Triton, CUDA, etc.).
For each subfolder:
1. Reads the query sequence and MSA from `sequences.a3m`.
2. Writes a `boltz_config.yaml` input file with a structural template if a
   `templates_complete/` subfolder is present, without otherwise.
3. Skips folders where `bz/predictions/sequences/models/` is non-empty.
4. Runs 5 seeds sequentially (from a shared base seed drawn once before the loop),
   each in its own `bz/seed_<s>/` output subfolder.
5. Reorganizes all seed outputs with `organize_files_boltz` after all seeds complete.
"""
function run_boltz2(
    clusters_folder::String,
    SIF_IMAGE_PATH::String,
    CACHE_DIR_Boltz::String,
)
    if isempty(clusters_folder)
        throw(ErrorException("No cluster_* folders were found in $clusters_folder"))
    end
    #=
    input_dirs = sort(
        filter(
            d -> isdir(d) && startswith(basename(d), "cluster_"),
            readdir(clusters_folder; join=true)
        )
    )
    =#
    input_dirs = sort(filter(d -> isdir(d), readdir(clusters_folder; join = true)))
    println("📂 Number of folder found: ", length(input_dirs))
    run_name = split(basename(clusters_folder), "_")[1]
    seed = rand(10_000:99_999)
    for input_dir in input_dirs
        #Get the input informations
        cluster_name = basename(input_dir)
        out_json=joinpath(input_dir, "boltz_config.yaml")

        ref_msa_path = joinpath(input_dir, "sequences.a3m")
        ids, sequences = read_a3m(ref_msa_path)
        query_id = ids[1]
        chain_id = split(query_id, "_")[2]
        sequence = sequences[1]
        if isdir(joinpath(input_dir, "templates_complete"))
            template_name=glob("*.cif", joinpath(input_dir, "templates_complete"))[1]
            template_path =
                joinpath("/mnt/input/templates_complete", basename(template_name))
            split_name = String(first(splitext(basename(template_name))))
            real_name = [
                ids[i] for i in eachindex(ids) if startswith(ids[i], lowercase(split_name))
            ]
            chain_template=split(real_name[1], "_")
            if length(chain_template) > 1
                chain_template = String(uppercase(chain_template[2]))
            else
                chain_template = "A"
            end
            # Right input file
            config = Dict(
                "sequences" => [
                    Dict(
                        "protein" => Dict(
                            "id" => chain_id,
                            "sequence" => sequence,
                            "msa" => "/mnt/input/sequences.a3m",
                        ),
                    ),
                ],
                "templates" =>
                    [Dict("cif" => template_path, "template_id" => chain_template)],
            )

            open(out_json, "w") do io
                YAML.write(io, config)
            end

        else
            # Right input file without template
            # Right input file
            config = Dict(
                "sequences" => [
                    Dict(
                        "protein" => Dict(
                            "id" => chain_id,
                            "sequence" => sequence,
                            "msa" => "/mnt/input/sequences.a3m",
                        ),
                    ),
                ],
            )

            open(out_json, "w") do io
                YAML.write(io, config)
            end

        end

        #Create output folder
        output_dir=joinpath(input_dir, "bz")
        if isdir(output_dir)
            if isdir(joinpath(output_dir, "predictions", "sequences", "models"))
                out=glob("*", joinpath(output_dir, "predictions", "sequences", "models"))
                @show out
                if !isempty(out)
                    @show "Folder already process"
                    continue
                else
                    rm(output_dir; recursive = true, force = true)
                end
            else
                rm(output_dir; recursive = true, force = true)
            end
        end
        mkdir(output_dir)

        println("🚀Start $cluster_name ...")
        for s = seed:(seed+4)
            # Run AlphaFold3
            output_dir_seed=joinpath(output_dir, "seed_$s")
            isdir(output_dir_seed) || mkdir(output_dir_seed)
            cmd = `apptainer exec --nv --no-home --cleanenv \
                --bind $input_dir:/mnt/input \
                --bind $output_dir_seed:/mnt/output \
                --bind $CACHE_DIR_Boltz:/cache \
                --env NUMBA_CACHE_DIR=/cache/numba \
                --env HOME=/cache \
                --env XDG_CACHE_HOME=/cache \
                --env TRITON_CACHE_DIR=/cache/triton \
                --env CUDA_CACHE_PATH=/cache/cuda \
                --env TORCHINDUCTOR_CACHE_DIR=/cache/torchinductor \
                --env CUEQUIVARIANCE_CACHE_DIR=/cache/cuequivariance \
                $SIF_IMAGE_PATH \
                    boltz predict /mnt/input/boltz_config.yaml \
                --out_dir /mnt/output/ \
                --cache /cache \
                --seed $s \
                --diffusion_samples 5 \
                --override  
                `
            try
                run(cmd)
                @info "✅ Boltz2 run finished for $cluster_name"
                #organize output files

            catch e
                @warn e
            end
        end
        organize_files_boltz(output_dir)
    end
end
