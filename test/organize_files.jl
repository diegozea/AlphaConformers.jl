using TestItems

@testitem "safe_move moves matching files and ignores empty matches" begin
    mktempdir() do dir
        dst = joinpath(dir, "dst")
        mkdir(dst)
        write(joinpath(dir, "a.txt"), "alpha")
        write(joinpath(dir, "b.dat"), "beta")

        AlphaConformers.safe_move(joinpath(dir, "*.txt"), dst)
        AlphaConformers.safe_move(joinpath(dir, "missing*"), dst)

        @test !isfile(joinpath(dir, "a.txt"))
        @test read(joinpath(dst, "a.txt"), String) == "alpha"
        @test isfile(joinpath(dir, "b.dat"))
    end
end

@testitem "organize_files creates ColabFold prediction layout" begin
    mktempdir() do dir
        write(joinpath(dir, "config.json"), "{\"config\":true}")
        write(joinpath(dir, "log.txt"), "log")
        write(joinpath(dir, "cite.bibtex"), "citation")
        write(joinpath(dir, "sequences.a3m"), ">query\nACD\n")
        write(joinpath(dir, "sequences.done.txt"), "done")
        write(joinpath(dir, "model_1.pdb"), "MODEL")
        write(joinpath(dir, "plot.png"), "PNG")
        write(joinpath(dir, "scores.json"), "{\"ptm\":0.9}")

        cwd = pwd()
        organize_files(dir)

        @test pwd() == cwd
        @test read(joinpath(dir, "predictions", "config.json"), String) ==
              "{\"config\":true}"
        @test read(joinpath(dir, "predictions", "log.txt"), String) == "log"
        @test read(joinpath(dir, "predictions", "cite.bibtex"), String) == "citation"
        @test isfile(joinpath(dir, "predictions", "sequences", "sequences.a3m"))
        @test isfile(joinpath(dir, "predictions", "sequences", "sequences.done.txt"))
        @test read(
            joinpath(dir, "predictions", "sequences", "models", "model_1.pdb"),
            String,
        ) == "MODEL"
        @test read(
            joinpath(dir, "predictions", "sequences", "plots", "plot.png"),
            String,
        ) == "PNG"
        @test read(
            joinpath(dir, "predictions", "sequences", "scores", "scores.json"),
            String,
        ) == "{\"ptm\":0.9}"
    end
end

@testitem "triage_outputs delegates to current prediction filter" begin
    mktempdir() do dir
        query_struct = joinpath(dir, "query.pdb")
        expected = AlphaConformers.found_best_prediction(dir, query_struct, nothing, "af")
        actual = AlphaConformers.triage_outputs(dir, query_struct, nothing, "af")

        @test actual == expected
        @test isempty(actual[1])
        @test isempty(actual[2])
    end
end

@testitem "organize_files_af3 copies AlphaFold3 outputs and seed files" begin
    mktempdir() do dir
        run_dir = joinpath(dir, "query")
        mkdir(run_dir)
        write(joinpath(run_dir, "query_data.json"), "{\"name\":\"query\"}")
        write(joinpath(run_dir, "TERMS_OF_USE.md"), "terms")
        write(joinpath(run_dir, "ranking_scores.csv"), "rank,score\n1,0.9\n")
        write(joinpath(run_dir, "query_model.cif"), "MODEL")
        write(joinpath(run_dir, "query_confidences.json"), "{\"c\":1}")
        write(joinpath(run_dir, "query_summary_confidences.json"), "{\"s\":1}")
        seed_dir = joinpath(run_dir, "seed_10")
        mkdir(seed_dir)
        write(joinpath(seed_dir, "confidences.json"), "{\"seed\":10}")
        write(joinpath(seed_dir, "summary_confidences.json"), "{\"summary\":10}")
        write(joinpath(seed_dir, "model.cif"), "SEEDMODEL")

        cwd = pwd()
        organize_files_af3(dir, "Query")

        @test pwd() == cwd
        @test read(joinpath(dir, "predictions", "config.json"), String) ==
              "{\"name\":\"query\"}"
        @test read(joinpath(dir, "predictions", "cite.bibtex"), String) == "terms"
        @test isfile(joinpath(dir, "predictions", "ranking_scores.csv"))
        @test read(joinpath(dir, "predictions", "sequences", "query_model.cif"), String) ==
              "MODEL"
        @test read(
            joinpath(dir, "predictions", "sequences", "scores", "query_confidences.json"),
            String,
        ) == "{\"c\":1}"
        @test read(
            joinpath(dir, "predictions", "sequences", "scores", "seed_10_confidences.json"),
            String,
        ) == "{\"seed\":10}"
        @test read(
            joinpath(dir, "predictions", "sequences", "models", "seed_10.cif"),
            String,
        ) == "SEEDMODEL"
    end
end

@testitem "organize_files_boltz copies seed outputs with seed prefixes" begin
    mktempdir() do dir
        seed_dir = joinpath(dir, "seed_123")
        output_dir = joinpath(seed_dir, "run_1")
        mkpath(joinpath(output_dir, "processed"))
        mkpath(joinpath(output_dir, "lightning_logs", "version_0"))
        config_dir = joinpath(output_dir, "predictions", "query")
        mkpath(config_dir)
        write(joinpath(output_dir, "processed", "manifest.json"), "{\"manifest\":true}")
        write(joinpath(output_dir, "lightning_logs", "version_0", "hparams.yaml"), "x: 1\n")
        write(joinpath(config_dir, "model.cif"), "MODEL")
        write(joinpath(config_dir, "confidence.json"), "{\"score\":1}")
        write(joinpath(config_dir, "plot.npz"), "PLOT")

        cwd = pwd()
        organize_files_boltz(dir)

        @test pwd() == cwd
        @test read(joinpath(dir, "predictions", "config.json"), String) ==
              "{\"manifest\":true}"
        @test read(joinpath(dir, "predictions", "hparams.yaml"), String) == "x: 1\n"
        @test read(
            joinpath(dir, "predictions", "sequences", "models", "seed_123_model.cif"),
            String,
        ) == "MODEL"
        @test read(
            joinpath(dir, "predictions", "sequences", "scores", "seed_123_confidence.json"),
            String,
        ) == "{\"score\":1}"
        @test read(
            joinpath(dir, "predictions", "sequences", "plots", "seed_123_plot.npz"),
            String,
        ) == "PLOT"
    end
end
