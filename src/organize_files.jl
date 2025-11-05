
# Crée un dossier s'il n'existe pas
function safe_mkdir(path::AbstractString)
    if !isdir(path)
        mkpath(path)
    end
end

# Déplace les fichiers correspondant à un motif vers un dossier
function safe_move(src_pattern::AbstractString, dst_dir::AbstractString)
    for src in glob(src_pattern)
        try
            mv(src, dst_dir; force=true)
            println("Moved: $src → $dst_dir")
        catch e
            println("⚠️ Could not move $src: $e")
        end
    end
end

# Fonction principale
function organize_files(output_dir::AbstractString)
    println("📂 Réorganisation des fichiers dans $output_dir ...")
    cd(output_dir)

    # Dossier principal
    predictions_dir = joinpath(output_dir, "predictions")
    mkdir(predictions_dir)

    mv("config.json", joinpath(predictions_dir, "config.json"))
    mv("log.txt", joinpath(predictions_dir, "log.txt"))
    mv("cite.bibtex", joinpath(predictions_dir, "cite.bibtex"))

    # Dossier sequences
    seq_dir = joinpath(predictions_dir, "sequences")
    mkdir(seq_dir)

    mv("sequences.a3m", joinpath(seq_dir, "sequences.a3m"))
    mv("sequences.done.txt", joinpath(seq_dir, "sequences.done.txt"))

    # Dossier models
    models_dir = joinpath(seq_dir, "models")
    mkdir(models_dir)

    for file in glob("*.pdb", output_dir)
        mv(file, joinpath(models_dir, basename(file)); force=true)
    end

    # Dossier plots
    plots_dir = joinpath(seq_dir, "plots")
    mkdir(plots_dir)
    for file in glob("*.png", output_dir)
        mv(file, joinpath(plots_dir, basename(file)); force=true)
    end

    # Dossier scores
    scores_dir = joinpath(seq_dir, "scores")
    mkdir(scores_dir)
    for file in glob("*.json", output_dir)
        mv(file, joinpath(scores_dir, basename(file)); force=true)
    end

    println("\n✅ Réorganisation terminée !")
end

