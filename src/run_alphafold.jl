"""
    run_alphafold(clusters_folder::String; colabfold_path::String=get(ENV, "COLABFOLD_PATH", ""))

This function runs AlphaFold 2 (colabfold_batch) for each cluster in the 
`clusters_folder`. The `colabfold_path` argument is the path to the colabfold_batch program.
Alternatively, you can set the COLABFOLD_PATH environment variable.
"""
function run_alphafold(clusters_folder::String; colabfold_path::String=get(ENV, "COLABFOLD_PATH", ""))
    if isempty(colabfold_path)
        throw(ErrorException("The path to ColabFold is not defined, please set the COLABFOLD_PATH environment variable or the colabfold_path keyword argument."))
    end
    cluster_folders = filter!(
        dir -> occursin("cluster_", dir), 
        readdir(clusters_folder, join=true))
    if isempty(cluster_folders)
        throw(ErrorException("No cluster_* folders were found in $clusters_folder"))
    end
    # remember the current working directory
    try
        # run AlphaFold for each cluster
        for folder in cluster_folders
            cd(folder)
            @info "Running AlphaFold for $(abspath(folder))"
            if isdir("templates")
                if !isempty(readdir("templates"))
                    af_command = `python3 $colabfold_path sequences.a3m af --use-templates 1 --msa-input --custom-template-path templates/ --num-seeds 5 --use-dropout --num-models 2 --overwrite-existing-results`
                    @info "Running AlphaFold command: $af_command"
                    success(run(af_command))
                else
                    @warn "No templates found in $(abspath("templates"))"
                end
            else
                @warn "There is no templates folder in $(abspath(folder))"
            end
        end
    finally
        # return to the original working directory
        cd(clusters_folder)
    end
end

# === FONCTIONS UTILITAIRES ===

function run_cmd(cmd::Cmd)
    try
        run(cmd)
        println("✅ end without error\n")
    catch e
        @warn "⚠️ Error of execution for : $e"
    end
end

function run_alphafold_one_run(clusters_folder::String, SIF_PATH::String, CACHE_DIR::String)
    ### Check input directories
    if isempty(clusters_folder)
        throw(ErrorException("No cluster_* folders were found in $clusters_folder"))
    end
    input_dirs = sort(
        filter(
            d -> isdir(d) && startswith(basename(d), "cluster_"),
            readdir(clusters_folder; join=true)
        )
    )

    println("📂 Number of folder found: ", length(input_dirs))
    
    for input_dir in input_dirs
        name = basename(input_dir)
        output_dir = joinpath(input_dir, "af")
        @show "Output directory: $output_dir"
        if isdir(output_dir)
            rm(output_dir; recursive=true, force=true)
        end
        mkdir(output_dir)

        println("🚀Start $name ...")

        cmd = `apptainer exec --nv --no-home --cleanenv \
            --bind $input_dir:/mnt/input \
            --bind $output_dir:/mnt/output \
            --bind $CACHE_DIR:/cache \
            $SIF_PATH \
            colabfold_batch /mnt/input/sequences.a3m /mnt/output \
                --custom-template-path /mnt/input/templates/ \
                --num-seeds 5 --use-dropout --num-models 2 --overwrite-existing-results`

        run_cmd(cmd)
         
        organize_files(output_dir)
    end

    println("🎉 All the ColabFold run are finish !")
end