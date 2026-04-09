
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
    if isfile("sequences.a3m") || isfile("sequences.done.txt")
        mv("sequences.a3m", joinpath(seq_dir, "sequences.a3m"))
        mv("sequences.done.txt", joinpath(seq_dir, "sequences.done.txt"))
    end
    
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

function organize_files_af3(output_dir::AbstractString,run_name::String)
    println("📂 Réorganisation des fichiers dans $output_dir ...")
    cd(output_dir)
    output_path=joinpath(output_dir,lowercase(run_name))
    # Dossier principal
    predictions_dir = joinpath(output_dir, "predictions")
    safe_mkdir(predictions_dir)
    cp(joinpath(output_path,lowercase(run_name)*"_data.json"), joinpath(predictions_dir, "config.json"); force=true)
    cp(joinpath(output_path,"TERMS_OF_USE.md"), joinpath(predictions_dir, "cite.bibtex"); force=true)
    cp(joinpath(output_path,"ranking_scores.csv"), joinpath(predictions_dir, "ranking_scores.csv"); force=true)

    seq_dir = joinpath(predictions_dir, "sequences")
    safe_mkdir(seq_dir)
    cp(joinpath(output_path,lowercase(run_name)*"_model.cif"), joinpath(seq_dir, lowercase(run_name)*"_model.cif"); force=true)

    scores_dir = joinpath(seq_dir, "scores")
    safe_mkdir(scores_dir)
    cp(joinpath(output_path,lowercase(run_name)*"_confidences.json"), joinpath(scores_dir, lowercase(run_name)*"_confidences.json"); force=true)
    cp(joinpath(output_path,lowercase(run_name)*"_summary_confidences.json"), joinpath(scores_dir, lowercase(run_name)*"_summary_confidences.json"); force=true)

    models_dir = joinpath(seq_dir, "models")
    safe_mkdir(models_dir)

    dirs=glob("seed*",output_path)
    for dir in dirs
        cp(joinpath(dir,"confidences.json"), joinpath(scores_dir,basename(dir)*"_confidences.json"); force=true)
        cp(joinpath(dir,"summary_confidences.json"), joinpath(scores_dir,basename(dir)*"_summary_confidences.json"); force=true)
        cp(joinpath(dir,"model.cif"), joinpath(models_dir,basename(dir)*".cif"); force=true)
    end

    println("\n✅ Réorganisation terminée !")
end

function organize_files_boltz(output_dir::AbstractString)
    println("📂 Réorganisation des fichiers dans $output_dir ...")
    cd(output_dir)
    
    # Dossier principal
    predictions_dir = joinpath(output_dir, "predictions")
    safe_mkdir(predictions_dir)

    seq_dir = joinpath(predictions_dir, "sequences")
    safe_mkdir(seq_dir)

    scores_dir = joinpath(seq_dir, "scores")
    safe_mkdir(scores_dir)

    models_dir = joinpath(seq_dir, "models")
    safe_mkdir(models_dir)

    plots_dir = joinpath(seq_dir, "plots")
    safe_mkdir(plots_dir)

    seed_dir=glob("seed*",output_dir)
    for seed in seed_dir
        output_path=glob("*",seed)[1]
        @show output_path
        cp(joinpath(output_path,"processed/manifest.json"), joinpath(predictions_dir, "config.json"); force=true)
        cp(joinpath(output_path,"lightning_logs/version_0/hparams.yaml"), joinpath(predictions_dir, "hparams.yaml"); force=true)

        config=glob("*",joinpath(output_path,"predictions"))[1]
        models=glob("*cif",config)
        
        for model in models
            cp(model, joinpath(models_dir,basename(seed)*"_"*basename(model)); force=true)
        end

        confidences=glob("*json",config)
        for confidence in confidences
            cp(confidence, joinpath(scores_dir,basename(seed)*"_"*basename(confidence)); force=true)
        end

        plots=glob("*npz",config)
        for plot in plots
            cp(plot, joinpath(plots_dir,basename(seed)*"_"*basename(plot)); force=true)
        end
    end

    println("\n✅ Réorganisation terminée !")
end