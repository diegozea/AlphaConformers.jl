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
    current_dir = pwd()
    try
        # run AlphaFold for each cluster
        for folder in cluster_folders
            cd(folder)
            @info "Running AlphaFold for $(abspath(folder))"
            if isdir("templates")
                if !isempty(readdir("templates"))
                    af_command = `$colabfold_path sequences.a3m af --use-templates 1 --msa-input --custom-template-path templates/ --num-seeds 5 --use-dropout --num-models 2 --overwrite-existing-results`
                    @info "Running AlphaFold command: $af_command"
                    run(af_command)
                else
                    @warn "No templates found in $(abspath("templates"))"
                end
            else
                @warn "There is no templates folder in $(abspath(folder))"
            end
        end
    finally
        # return to the original working directory
        cd(current_dir)
    end
end