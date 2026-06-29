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

@testitem "triage_outputs reports missing Foldseek results clearly" begin
    mktempdir() do dir
        query_struct = joinpath(dir, "query.pdb")

        err = try
            AlphaConformers.triage_outputs(dir, query_struct, nothing; folder_af2_result = "af")
        catch err
            err
        end

        @test err isa ArgumentError
        @test occursin("Foldseek", sprint(showerror, err))
        @test occursin("foldseek_results_folder", sprint(showerror, err))
    end
end

@testitem "triage Foldseek result folder discovery is explicit" begin
    function write_m8(path)
        row = [
            "query",
            "target",
            "1.0",
            "10",
            "0",
            "0",
            "1",
            "10",
            "1",
            "10",
            "1e-10",
            "100",
            "0.9",
            "0.9",
            "0.9",
            "0.1",
            "1.0",
        ]
        write(path, join(row, '\t') * "\n")
    end

    mktempdir() do dir
        output_dir = joinpath(dir, "run")
        result_dir = joinpath(output_dir, "custom_fullpdb_results")
        mkpath(result_dir)
        write_m8(joinpath(result_dir, "query_results.m8"))

        @test AlphaConformers._resolve_foldseek_results_folder(output_dir) ==
              abspath(result_dir)
        @test AlphaConformers._resolve_foldseek_results_folder(
            output_dir;
            foldseek_results_folder = "custom_fullpdb_results",
        ) == abspath(result_dir)

        cd(dir) do
            relative_output_dir = joinpath("outputs", "1ABC_A")
            documented_result_dir = joinpath(relative_output_dir, "fullpdb_results")
            mkpath(documented_result_dir)
            write_m8(joinpath(documented_result_dir, "query_results.m8"))
            mkpath("fullpdb_results")
            write_m8(joinpath("fullpdb_results", "query_results.m8"))

            @test AlphaConformers._resolve_foldseek_results_folder(
                relative_output_dir;
                foldseek_results_folder = joinpath(
                    relative_output_dir,
                    "fullpdb_results",
                ),
            ) == abspath(documented_result_dir)

            @test AlphaConformers._resolve_foldseek_results_folder(
                relative_output_dir;
                foldseek_results_folder = "fullpdb_results",
            ) == abspath(documented_result_dir)
        end

        empty_dir = joinpath(output_dir, "empty_results")
        mkpath(empty_dir)
        err = try
            AlphaConformers._resolve_foldseek_results_folder(
                output_dir;
                foldseek_results_folder = empty_dir,
            )
        catch err
            err
        end

        @test err isa ArgumentError
        @test occursin("does not contain any `.m8`", sprint(showerror, err))

        second_dir = joinpath(output_dir, "other_results")
        mkpath(second_dir)
        write_m8(joinpath(second_dir, "query_results.m8"))

        err = try
            AlphaConformers._resolve_foldseek_results_folder(output_dir)
        catch err
            err
        end

        @test err isa ArgumentError
        @test occursin("multiple Foldseek result folders", sprint(showerror, err))
        @test occursin("foldseek_results_folder", sprint(showerror, err))
    end
end

@testitem "triage search results include rows from target_db_results" begin
    import DataFrames

    full_search_results = DataFrames.DataFrame(
        query = ["query"],
        target = ["1ABC_A"],
        qstart = [1],
        qend = [10],
        tstart = [1],
        tend = [10],
    )
    target_search_results = DataFrames.DataFrame(
        query = ["query"],
        target = ["2DEF_A"],
        qstart = [1],
        qend = [10],
        tstart = [1],
        tend = [10],
    )

    merged = AlphaConformers._merge_triage_search_results(
        full_search_results,
        target_search_results,
    )

    @test sort(merged.target) == ["1ABC_A", "2DEF_A"]
end

@testitem "triage UniProt lookup uses the Foldseek target chain" begin
    import DataFrames

    sifts_mapping = DataFrames.DataFrame(
        PDB = ["4lyl", "4lyl", "1akz", "6vba"],
        CHAIN = ["G", "H", "A", "A"],
        SP_PRIMARY = ["Q9I983", "P14739", "Q9I983", "Q9I983"],
    )
    search_results = DataFrames.DataFrame(
        query = fill("1AKZ_A", 12),
        target = [
            "4LYL_G",
            "4LYL_H",
            "4LYL.cif.gz_H",
            "4LYL_H.pdb",
            "1AKZ_A",
            "1AKZ_A.pdb",
            "1AKZ.cif",
            "1AKZ.cif.cif_A",
            "6VBA_A",
            "6VBA_A.cif.gz",
            "6VBA.cif",
            "6VBA.cif.gz_A",
        ],
    )

    result = AlphaConformers.found_uniprot_structure(
        search_results,
        sifts_mapping,
        "4LYL",
        "G",
    )

    @test Set(result.target) == Set([
        "4LYL_G",
        "1AKZ_A",
        "1AKZ_A.pdb",
        "1AKZ.cif",
        "1AKZ.cif.cif_A",
        "6VBA_A",
        "6VBA_A.cif.gz",
        "6VBA.cif",
        "6VBA.cif.gz_A",
    ])
    @test !("4LYL_H" in result.target)
    @test !("4LYL.cif.gz_H" in result.target)
    @test !("4LYL_H.pdb" in result.target)
end

@testitem "triage UniProt lookup normalizes query Foldseek names" begin
    import DataFrames

    sifts_mapping = DataFrames.DataFrame(
        PDB = ["4lyl", "1akz", "6vba"],
        CHAIN = ["G", "A", "A"],
        SP_PRIMARY = ["Q9I983", "Q9I983", "Q9I983"],
    )
    search_results = DataFrames.DataFrame(
        query = fill("1AKZ_A", 5),
        target = ["4LYL_G", "1AKZ_A.pdb", "1AKZ.cif", "6VBA_A.cif.gz", "6VBA.cif"],
    )
    expected_targets =
        Set(["4LYL_G", "1AKZ_A.pdb", "1AKZ.cif", "6VBA_A.cif.gz", "6VBA.cif"])

    chain_suffix_result = AlphaConformers.found_uniprot_structure(
        search_results,
        sifts_mapping,
        "1AKZ",
        "A.pdb",
    )
    @test Set(chain_suffix_result.target) == expected_targets

    pdb_and_chain_suffix_result = AlphaConformers.found_uniprot_structure(
        search_results,
        sifts_mapping,
        "6VBA.cif",
        "A.cif.gz",
    )
    @test Set(pdb_and_chain_suffix_result.target) == expected_targets
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
