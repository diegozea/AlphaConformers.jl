using TestItems

@testitem "Pipeline helper functions parse names, progress, and sequence projections" begin
    import MIToS
    import Printf

    function write_numbered_pdb(path)
        open(path, "w") do io
            for (serial, resnum) in enumerate((1, 2, 10))
                Printf.@printf(
                    io,
                    "ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n",
                    serial,
                    "CA",
                    "ALA",
                    resnum,
                    1.0 * serial,
                    0.0,
                    0.0,
                    "C",
                )
            end
            println(io, "END")
        end
    end

    @test AlphaConformers._parse_query_pdb_filename("/tmp/1ABC_A.pdb") == ("1ABC", "A")
    @test AlphaConformers._parse_query_pdb_filename("AF-P12345-F1-model_v4_A.pdb") ==
          ("AF-P12345-F1-model", "v4")
    @test_throws ErrorException AlphaConformers._parse_query_pdb_filename("1ABC.pdb")
    @test_throws ErrorException AlphaConformers._parse_query_pdb_filename("_A.pdb")

    @test AlphaConformers._aligned_pdb_file("query", "target", "/tmp/aligned") ==
          joinpath("/tmp/aligned", "aln_query_target.pdb")

    progress = AlphaConformers.Progress(["one", "two"], 0, 2)
    redirect_stdout(devnull) do
        AlphaConformers.progress_bar(progress, "one")
    end
    @test progress.current == 1

    steps = [
        "Foldseek search",
        "Merging tables",
        "Adding known conformations",
        "Merging MSAs",
        "Final MSA",
        "Aligned structures",
        "Clustering Hobohm",
        "Creating folders + cleanup",
    ]
    final_progress = AlphaConformers.Progress(steps, 0, length(steps))
    redirect_stdout(devnull) do
        for step in steps
            AlphaConformers.progress_bar(final_progress, step)
        end
    end
    @test final_progress.current == length(steps)

    projected = AlphaConformers.align_full_seq("XABCY", ["ABC", "DEF", "D-F"])
    @test projected == ["XABCY", "-DEF-", "-D-F-"]

    msa = AlphaConformers.vector_to_msa(["ACD", "A-D"])
    @test size(msa) == (2, 3)
    @test MIToS.MSA.sequencenames(msa) == ["seq_1", "seq_2"]

    prepared = AlphaConformers.PreparedInputs("query.pdb", "run", "run/fullpdb_results")
    @test prepared.query_struct == "query.pdb"
    @test prepared.output_dir == "run"
    @test prepared.foldseek_results_folder == "run/fullpdb_results"

    mktempdir() do dir
        pdb = joinpath(dir, "numbered.pdb")
        write_numbered_pdb(pdb)
        @test AlphaConformers.get_res_beg_end(pdb) == (1, 10)
    end
end

@testitem "alphaconformers runs valid step combinations and forwards predictor arguments" begin
    import DataFrames

    function keyword_dict(kwargs)
        dict = Dict{Symbol,Any}()
        for pair in kwargs
            dict[pair.first] = pair.second
        end
        return dict
    end

    function fake_steps(calls)
        fake_preparer = function (query_struct, pdb_folder, output_dir; kwargs...)
            calls[:prepare] = (
                query_struct = query_struct,
                pdb_folder = pdb_folder,
                output_dir = output_dir,
                kwargs = keyword_dict(kwargs),
            )
            return :prepared
        end
        fake_predictor = function (output_dir, args...; kwargs...)
            calls[:predict] =
                (output_dir = output_dir, args = args, kwargs = keyword_dict(kwargs))
            return :predicted
        end
        fake_triager = function (output_dir, query_struct, sifts_uniprot_mapping; kwargs...)
            calls[:triage] = (
                output_dir = output_dir,
                query_struct = query_struct,
                sifts_uniprot_mapping = sifts_uniprot_mapping,
                kwargs = keyword_dict(kwargs),
            )
            return (:triaged, output_dir)
        end
        fake_deepaccnet = function (output_dir, sif_path; kwargs...)
            calls[:deepaccnet] = (
                output_dir = output_dir,
                sif_path = sif_path,
                kwargs = keyword_dict(kwargs),
            )
            csv_path = joinpath(output_dir, "deepaccnet_results.csv")
            write(csv_path, "sample,cb-lddt\ncluster_1_m1,0.9\ncluster_1_m2,0.6\n")
            return csv_path
        end
        fake_cluster = function (output_dir, ref1, ref2; kwargs...)
            calls[:cluster] = (
                output_dir = output_dir,
                ref1 = ref1,
                ref2 = ref2,
                kwargs = keyword_dict(kwargs),
            )
            return :clustered
        end
        return fake_preparer, fake_predictor, fake_triager, fake_deepaccnet, fake_cluster
    end

    combinations = [
        (name = "prepare only", prepare = true, predict = false, triage = false),
        (name = "prepare and predict", prepare = true, predict = true, triage = false),
        (name = "full pipeline", prepare = true, predict = true, triage = true),
        (name = "predict only", prepare = false, predict = true, triage = false),
        (name = "predict and triage", prepare = false, predict = true, triage = true),
        (name = "triage only", prepare = false, predict = false, triage = true),
    ]

    mktempdir() do dir
        query_struct = joinpath(dir, "1ABC_A.pdb")
        pdb_folder = joinpath(dir, "pdb")
        output_dir = joinpath(dir, "run")
        mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[])

        models = joinpath(
            output_dir,
            "cluster_1",
            "custom_af",
            "predictions",
            "sequences",
            "models",
        )
        mkpath(models)
        write(joinpath(models, "m1.pdb"), "")
        write(joinpath(models, "m2.pdb"), "")

        for combination in combinations
            calls = Dict{Symbol,Any}()
            fake_preparer, fake_predictor, fake_triager, fake_deepaccnet, fake_cluster =
                fake_steps(calls)

            result = redirect_stdout(devnull) do
                if combination.predict
                    AlphaConformers._alphaconformers(
                        fake_preparer,
                        fake_triager,
                        fake_deepaccnet,
                        fake_cluster,
                        "image.sif",
                        "cache";
                        query_struct,
                        pdb_folder,
                        output_dir,
                        prepare = combination.prepare,
                        predict = combination.predict,
                        triage = combination.triage,
                        structure_predictor = fake_predictor,
                        databases = ["foldseek_db"],
                        sifts_uniprot_mapping = mapping,
                        folder_af2_result = "custom_af",
                        deepaccnet_sif = "image.sif",
                        predictor_keyword = 42,
                    )
                else
                    AlphaConformers._alphaconformers(
                        fake_preparer,
                        fake_triager,
                        fake_deepaccnet,
                        fake_cluster;
                        query_struct,
                        pdb_folder,
                        output_dir,
                        prepare = combination.prepare,
                        predict = combination.predict,
                        triage = combination.triage,
                        structure_predictor = fake_predictor,
                        databases = ["foldseek_db"],
                        sifts_uniprot_mapping = mapping,
                        folder_af2_result = "custom_af",
                        deepaccnet_sif = "image.sif",
                    )
                end
            end

            @test result.output_dir == abspath(output_dir)
            @test result.prepared == combination.prepare
            @test result.predicted == combination.predict
            @test result.triaged == combination.triage
            @test (result.triage_result !== nothing) == combination.triage
            # DeepAccNet and clustering are sub-steps of triage: they run exactly when it does.
            @test (result.deepaccnet_result !== nothing) == combination.triage
            @test (result.cluster_result !== nothing) == combination.triage
            @test haskey(calls, :prepare) == combination.prepare
            @test haskey(calls, :predict) == combination.predict
            @test haskey(calls, :triage) == combination.triage
            @test haskey(calls, :deepaccnet) == combination.triage
            @test haskey(calls, :cluster) == combination.triage

            if combination.prepare
                @test calls[:prepare].query_struct == abspath(query_struct)
                @test calls[:prepare].pdb_folder == abspath(pdb_folder)
                @test calls[:prepare].output_dir == abspath(output_dir)
                @test calls[:prepare].kwargs[:databases] == ["foldseek_db"]
                @test calls[:prepare].kwargs[:sifts_uniprot_mapping] === mapping
                @test !haskey(calls[:prepare].kwargs, :predictor_keyword)
            end

            if combination.predict
                @test calls[:predict].output_dir == abspath(output_dir)
                @test calls[:predict].args == ("image.sif", "cache")
                @test calls[:predict].kwargs[:predictor_keyword] == 42
            end

            if combination.triage
                @test calls[:triage].output_dir == abspath(output_dir)
                @test calls[:triage].query_struct == abspath(query_struct)
                @test calls[:triage].sifts_uniprot_mapping === mapping
                @test calls[:triage].kwargs[:folder_af2_result] == "custom_af"
                @test !haskey(calls[:triage].kwargs, :predictor_keyword)

                @test calls[:deepaccnet].output_dir == abspath(output_dir)
                @test calls[:deepaccnet].sif_path == "image.sif"
                @test calls[:deepaccnet].kwargs[:n_threads] isa Int

                @test calls[:cluster].output_dir == abspath(output_dir)
                @test calls[:cluster].kwargs[:method] == "custom_af"
            end
        end
    end
end

@testitem "alphaconformers public wrapper forwards predictor arguments only to prediction" begin
    function keyword_dict(kwargs)
        dict = Dict{Symbol,Any}()
        for pair in kwargs
            dict[pair.first] = pair.second
        end
        return dict
    end

    calls = Ref{Any}()
    fake_predictor = function (output_dir, args...; kwargs...)
        calls[] = (output_dir = output_dir, args = args, kwargs = keyword_dict(kwargs))
        return :predicted
    end

    mktempdir() do dir
        result = redirect_stdout(devnull) do
            AlphaConformers.alphaconformers(
                "image.sif",
                "cache";
                output_dir = dir,
                prepare = false,
                predict = true,
                triage = false,
                structure_predictor = fake_predictor,
                predictor_keyword = 42,
            )
        end

        @test result.output_dir == abspath(dir)
        @test result.prepared == false
        @test result.predicted == true
        @test result.triage_result === nothing
        @test calls[].output_dir == abspath(dir)
        @test calls[].args == ("image.sif", "cache")
        @test calls[].kwargs[:predictor_keyword] == 42
    end
end

@testitem "alphaconformers passes prepared Foldseek result folder to triage" begin
    import DataFrames

    function keyword_dict(kwargs)
        dict = Dict{Symbol,Any}()
        for pair in kwargs
            dict[pair.first] = pair.second
        end
        return dict
    end

    calls = Dict{Symbol,Any}()
    fake_preparer = function (query_struct, pdb_folder, output_dir; kwargs...)
        folder = joinpath(output_dir, "custom_fullpdb_results")
        calls[:prepare] = (
            query_struct = query_struct,
            pdb_folder = pdb_folder,
            output_dir = output_dir,
            kwargs = keyword_dict(kwargs),
        )
        return (;
            query_struct = query_struct,
            output_dir = output_dir,
            foldseek_results_folder = folder,
        )
    end
    fake_predictor = function (output_dir, args...; kwargs...)
        calls[:predict] =
            (output_dir = output_dir, args = args, kwargs = keyword_dict(kwargs))
        return :predicted
    end
    fake_triager = function (
        output_dir,
        query_struct,
        sifts_uniprot_mapping;
        foldseek_results_folder,
        kwargs...,
    )
        calls[:triage] = (
            output_dir = output_dir,
            query_struct = query_struct,
            sifts_uniprot_mapping = sifts_uniprot_mapping,
            foldseek_results_folder = foldseek_results_folder,
            kwargs = keyword_dict(kwargs),
        )
        return :triaged
    end
    fake_deepaccnet = function (output_dir, sif_path; kwargs...)
        csv_path = joinpath(output_dir, "deepaccnet_results.csv")
        write(csv_path, "sample,cb-lddt\ncluster_1_m1,0.9\ncluster_1_m2,0.6\n")
        return csv_path
    end
    fake_cluster = (output_dir, ref1, ref2; kwargs...) -> :clustered

    mktempdir() do dir
        query_struct = joinpath(dir, "1ABC_A.pdb")
        pdb_folder = joinpath(dir, "pdb")
        output_dir = joinpath(dir, "run")
        mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[])

        # Minimal prediction tree so the triage step's DeepAccNet-to-clustering adapter resolves.
        models =
            joinpath(output_dir, "cluster_1", "af", "predictions", "sequences", "models")
        mkpath(models)
        write(joinpath(models, "m1.pdb"), "")
        write(joinpath(models, "m2.pdb"), "")

        result = redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                fake_preparer,
                fake_triager,
                fake_deepaccnet,
                fake_cluster,
                "image.sif",
                "cache";
                query_struct,
                pdb_folder,
                output_dir,
                prepare = true,
                predict = true,
                triage = true,
                structure_predictor = fake_predictor,
                databases = [joinpath(dir, "custom_fullpdb")],
                sifts_uniprot_mapping = mapping,
                deepaccnet_sif = "image.sif",
            )
        end

        expected_folder = joinpath(abspath(output_dir), "custom_fullpdb_results")
        @test result.triage_result == :triaged
        @test calls[:triage].foldseek_results_folder == expected_folder

        explicit_folder = joinpath(abspath(output_dir), "chosen_results")
        result = redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                fake_preparer,
                fake_triager,
                fake_deepaccnet,
                fake_cluster,
                "image.sif",
                "cache";
                query_struct,
                pdb_folder,
                output_dir,
                prepare = true,
                predict = true,
                triage = true,
                structure_predictor = fake_predictor,
                databases = [joinpath(dir, "custom_fullpdb")],
                sifts_uniprot_mapping = mapping,
                foldseek_results_folder = explicit_folder,
                deepaccnet_sif = "image.sif",
            )
        end

        @test result.triage_result == :triaged
        @test calls[:triage].foldseek_results_folder == explicit_folder
    end
end

@testitem "alphaconformers selects the PDB Foldseek result folder before preparation" begin
    import DataFrames

    calls = Dict{Symbol,Any}()
    fake_preparer =
        function (query_struct, pdb_folder, output_dir; foldseek_results_folder, kwargs...)
            calls[:prepare] = (foldseek_results_folder = foldseek_results_folder,)
            return AlphaConformers.PreparedInputs(
                query_struct,
                output_dir,
                foldseek_results_folder,
            )
        end
    fake_predictor = function (output_dir, args...; kwargs...)
        calls[:predict] = true
        return :predicted
    end
    fake_triager = function (
        output_dir,
        query_struct,
        sifts_uniprot_mapping;
        foldseek_results_folder,
        kwargs...,
    )
        calls[:triage] = (foldseek_results_folder = foldseek_results_folder,)
        return :triaged
    end
    fake_deepaccnet = function (output_dir, sif_path; kwargs...)
        csv_path = joinpath(output_dir, "deepaccnet_results.csv")
        write(csv_path, "sample,cb-lddt\ncluster_1_m1,0.9\ncluster_1_m2,0.6\n")
        return csv_path
    end
    fake_cluster = (output_dir, ref1, ref2; kwargs...) -> :clustered

    mktempdir() do dir
        query_struct = joinpath(dir, "1ABC_A.pdb")
        pdb_folder = joinpath(dir, "pdb")
        output_dir = joinpath(dir, "run")
        mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[])

        models =
            joinpath(output_dir, "cluster_1", "af", "predictions", "sequences", "models")
        mkpath(models)
        write(joinpath(models, "m1.pdb"), "")
        write(joinpath(models, "m2.pdb"), "")

        result = redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                fake_preparer,
                fake_triager,
                fake_deepaccnet,
                fake_cluster,
                "image.sif",
                "cache";
                query_struct,
                pdb_folder,
                output_dir,
                prepare = true,
                predict = true,
                triage = true,
                structure_predictor = fake_predictor,
                databases = [joinpath(dir, "afdb"), joinpath(dir, "myPDBdatabase")],
                sifts_uniprot_mapping = mapping,
                deepaccnet_sif = "image.sif",
            )
        end

        expected_folder = joinpath(abspath(output_dir), "myPDBdatabase_results")
        @test result.triage_result == :triaged
        @test calls[:prepare].foldseek_results_folder == expected_folder
        @test calls[:triage].foldseek_results_folder == expected_folder
    end
end

@testitem "alphaconformers reports invalid step and keyword combinations clearly" begin
    import DataFrames

    function argument_error_message(f)
        try
            f()
        catch err
            @test err isa ArgumentError
            return sprint(showerror, err)
        end
        error("Expected ArgumentError")
    end

    should_not_run = function (args...; kwargs...)
        error("This test should fail during validation before running a step.")
    end
    mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[])

    mktempdir() do dir
        query_struct = joinpath(dir, "1ABC_A.pdb")
        pdb_folder = joinpath(dir, "pdb")
        output_dir = joinpath(dir, "run")

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                prepare = false,
                predict = false,
                triage = true,
                query_struct,
                sifts_uniprot_mapping = mapping,
            )
        end
        @test occursin("`output_dir` is required", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                prepare = false,
                predict = false,
                triage = false,
            )
        end
        @test occursin("At least one pipeline step must be enabled", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                pdb_folder,
                databases = ["foldseek_db"],
                prepare = true,
                predict = false,
                triage = false,
            )
        end
        @test occursin("`query_struct` is required when `prepare=true`", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                sifts_uniprot_mapping = mapping,
                prepare = false,
                predict = false,
                triage = true,
            )
        end
        @test occursin("`query_struct` is required when `triage=true`", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                query_struct,
                databases = ["foldseek_db"],
                prepare = true,
                predict = false,
                triage = false,
            )
        end
        @test occursin("`pdb_folder` is required when `prepare=true`", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                query_struct,
                pdb_folder,
                prepare = true,
                predict = false,
                triage = false,
            )
        end
        @test occursin("`databases` is required when `prepare=true`", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                query_struct,
                pdb_folder,
                databases = ["foldseek_db"],
                sifts_uniprot_mapping = mapping,
                prepare = true,
                predict = false,
                triage = true,
            )
        end
        @test occursin("`triage=true` requires prediction outputs", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                query_struct,
                prepare = false,
                predict = false,
                triage = true,
            )
        end
        @test occursin("`sifts_uniprot_mapping` is required", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run,
                "image.sif";
                output_dir,
                query_struct,
                pdb_folder,
                databases = ["foldseek_db"],
                prepare = true,
                predict = false,
                triage = false,
            )
        end
        @test occursin("Predictor arguments were provided but `predict=false`", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                query_struct,
                pdb_folder,
                databases = ["foldseek_db"],
                prepare = true,
                predict = false,
                triage = false,
                predictor_keyword = 42,
            )
        end
        @test occursin("Predictor arguments were provided but `predict=false`", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                should_not_run,
                should_not_run,
                should_not_run,
                should_not_run;
                output_dir,
                query_struct,
                sifts_uniprot_mapping = mapping,
                prepare = false,
                predict = false,
                triage = true,
            )
        end
        @test occursin("deepaccnet_sif", msg)

        msg = argument_error_message() do
            AlphaConformers.prepare_inputs(query_struct, pdb_folder, output_dir)
        end
        @test occursin("`databases` is required when running `prepare_inputs`", msg)

        prepare_called = Ref(false)
        fail_if_prepare_runs = function (args...; kwargs...)
            prepare_called[] = true
            error("Preparation should not run when the PDB Foldseek database is unclear.")
        end
        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                fail_if_prepare_runs,
                should_not_run,
                should_not_run,
                should_not_run,
                "image.sif",
                "cache";
                output_dir,
                query_struct,
                pdb_folder,
                databases = [joinpath(dir, "afdb"), joinpath(dir, "uniref")],
                sifts_uniprot_mapping = mapping,
                prepare = true,
                predict = true,
                triage = true,
                structure_predictor = should_not_run,
                deepaccnet_sif = "image.sif",
            )
        end
        @test prepare_called[] == false
        @test occursin("foldseek_results_folder", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                fail_if_prepare_runs,
                should_not_run,
                should_not_run,
                should_not_run,
                "image.sif",
                "cache";
                output_dir,
                query_struct,
                pdb_folder,
                databases = [joinpath(dir, "pdb_parent", "afdb"), joinpath(dir, "uniref")],
                sifts_uniprot_mapping = mapping,
                prepare = true,
                predict = true,
                triage = true,
                structure_predictor = should_not_run,
                deepaccnet_sif = "image.sif",
            )
        end
        @test prepare_called[] == false
        @test occursin("foldseek_results_folder", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                fail_if_prepare_runs,
                should_not_run,
                should_not_run,
                should_not_run,
                "image.sif",
                "cache";
                output_dir,
                query_struct,
                pdb_folder,
                databases = [joinpath(dir, "fullpdb"), joinpath(dir, "pdb")],
                sifts_uniprot_mapping = mapping,
                prepare = true,
                predict = true,
                triage = true,
                structure_predictor = should_not_run,
                deepaccnet_sif = "image.sif",
            )
        end
        @test prepare_called[] == false
        @test occursin("foldseek_results_folder", msg)
    end
end

@testitem "alphaconformers runs the triage step's DeepAccNet and clustering sub-steps" begin
    import DataFrames

    function keyword_dict(kwargs)
        dict = Dict{Symbol,Any}()
        for pair in kwargs
            dict[pair.first] = pair.second
        end
        return dict
    end

    mktempdir() do dir
        output_dir = joinpath(dir, "run")
        # A minimal prediction tree so the real DeepAccNet-to-clustering adapter resolves.
        models =
            joinpath(output_dir, "cluster_1", "af", "predictions", "sequences", "models")
        mkpath(models)
        write(joinpath(models, "m1.pdb"), "")
        write(joinpath(models, "m2.pdb"), "")
        query_struct = joinpath(dir, "1ABC_A.pdb")
        mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[])

        calls = Dict{Symbol,Any}()
        fake_triager = function (output_dir, query_struct, sifts_uniprot_mapping; kwargs...)
            calls[:triage] = (output_dir = output_dir, kwargs = keyword_dict(kwargs))
            return :triaged
        end
        fake_deepaccnet = function (output_dir, sif_path; kwargs...)
            calls[:deepaccnet] = (
                output_dir = output_dir,
                sif_path = sif_path,
                kwargs = keyword_dict(kwargs),
            )
            csv_path = joinpath(output_dir, "deepaccnet_results.csv")
            write(csv_path, "sample,cb-lddt\ncluster_1_m1,0.9\ncluster_1_m2,0.6\n")
            return csv_path
        end
        fake_cluster = function (output_dir, ref1, ref2; kwargs...)
            calls[:cluster] = (
                output_dir = output_dir,
                ref1 = ref1,
                ref2 = ref2,
                kwargs = keyword_dict(kwargs),
            )
            return :clustered
        end

        # Analysis-only run: triage drives its three sub-steps in order against an existing
        # folder, feeding the DeepAccNet scores into clustering.
        result = redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                (args...; kwargs...) -> error("prepare must not run"),
                fake_triager,
                fake_deepaccnet,
                fake_cluster;
                output_dir,
                prepare = false,
                predict = false,
                triage = true,
                query_struct,
                sifts_uniprot_mapping = mapping,
                deepaccnet_sif = "image.sif",
                ref1 = "apo.pdb",
                ref2 = "holo.pdb",
            )
        end

        @test result.output_dir == abspath(output_dir)
        @test result.triaged == true
        @test result.triage_result == :triaged
        @test result.deepaccnet_result !== nothing
        @test result.cluster_result == :clustered
        @test haskey(calls, :triage)
        @test haskey(calls, :deepaccnet)
        @test haskey(calls, :cluster)

        @test calls[:deepaccnet].output_dir == abspath(output_dir)
        @test calls[:deepaccnet].sif_path == "image.sif"
        @test calls[:deepaccnet].kwargs[:n_threads] isa Int
        @test calls[:cluster].output_dir == abspath(output_dir)
        @test calls[:cluster].ref1 == "apo.pdb"
        @test calls[:cluster].ref2 == "holo.pdb"
        @test calls[:cluster].kwargs[:method] == "af"

        # With triage disabled, none of its three sub-steps run.
        no_run =
            (args...; kwargs...) -> error("triage sub-steps must not run when triage=false")
        predict_only = redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                (args...; kwargs...) -> error("prepare must not run"),
                no_run,
                no_run,
                no_run;
                output_dir,
                prepare = false,
                predict = true,
                triage = false,
                structure_predictor = (output_dir, args...; kwargs...) -> :predicted,
            )
        end
        @test predict_only.predicted == true
        @test predict_only.triaged == false
        @test predict_only.triage_result === nothing
        @test predict_only.deepaccnet_result === nothing
        @test predict_only.cluster_result === nothing
    end
end

@testitem "alphaconformers feeds DeepAccNet scores into clustering" begin
    import DataFrames

    no_prepare = (args...; kwargs...) -> error("prepare must not run")
    fake_triage = (args...; kwargs...) -> :triaged

    # A fake DeepAccNet that writes a real results CSV keyed by the flat symlink stem.
    fake_deepaccnet = function (output_dir, sif_path; kwargs...)
        csv_path = joinpath(output_dir, "deepaccnet_results.csv")
        write(csv_path, "sample,cb-lddt\ncluster_1_m1,0.9\ncluster_1_m2,0.6\n")
        return csv_path
    end

    mktempdir() do dir
        output_dir = joinpath(dir, "run")
        models =
            joinpath(output_dir, "cluster_1", "af", "predictions", "sequences", "models")
        mkpath(models)
        write(joinpath(models, "m1.pdb"), "")
        write(joinpath(models, "m2.pdb"), "")
        query_struct = joinpath(dir, "1ABC_A.pdb")
        mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[])

        # Auto-feed: with no explicit score table, the DeepAccNet CSV is converted and
        # handed to clustering as a relative-path-keyed DataFrame.
        captured = Dict{Symbol,Any}()
        fake_cluster = function (output_dir, ref1, ref2; score_table, kwargs...)
            captured[:score_table] = score_table
            return :clustered
        end

        redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                no_prepare,
                fake_triage,
                fake_deepaccnet,
                fake_cluster;
                output_dir,
                prepare = false,
                predict = false,
                triage = true,
                query_struct,
                sifts_uniprot_mapping = mapping,
                deepaccnet_sif = "image.sif",
            )
        end

        table = captured[:score_table]
        @test table isa DataFrames.DataFrame
        @test table.Name == [
            joinpath("cluster_1", "af", "predictions", "sequences", "models", "m1.pdb"),
            joinpath("cluster_1", "af", "predictions", "sequences", "models", "m2.pdb"),
        ]
        @test table[!, "cb-lddt"] == [0.9, 0.6]
    end
end

@testitem "alphaconformers runs real clustering fed by DeepAccNet scores" begin
    import Random
    import DataFrames

    ENV["GKSwstype"] = "100"  # headless GKS backend so figures render without a display

    # Same deterministic three-shape ensemble as the clustering fixtures.
    function _write_fixture(data_root, system)
        models = joinpath(
            data_root,
            system,
            "cluster_1",
            "af",
            "predictions",
            "sequences",
            "models",
        )
        mkpath(models)
        residues = 10
        base = [(3.8 * (i - 1), 0.0, 0.0) for i = 1:residues]
        for (group, bend) in ["A" => 0.0, "B" => 5.0, "C" => 12.0]
            for rep = 1:3
                rng = Random.MersenneTwister(abs(hash((group, rep))) % 100_000)
                path = joinpath(models, "conf_$(group)$(rep).pdb")
                open(path, "w") do io
                    for (i, b) in enumerate(base)
                        yb =
                            i > residues ÷ 2 ? bend * (i - residues ÷ 2) / (residues ÷ 2) :
                            0.0
                        jitter = 0.05 .* (Random.rand(rng, 3) .- 0.5)
                        x, y, z = b[1] + jitter[1], b[2] + yb + jitter[2], b[3] + jitter[3]
                        line = string(
                            "ATOM  ",
                            lpad(i, 5),
                            "  CA  ALA A",
                            lpad(i, 4),
                            "    ",
                            lpad(string(round(x, digits = 3)), 8),
                            lpad(string(round(y, digits = 3)), 8),
                            lpad(string(round(z, digits = 3)), 8),
                            "  1.00  0.00           C",
                        )
                        println(io, rpad(line, 80))
                    end
                    println(io, "END")
                end
            end
        end
        return models
    end

    mktempdir() do fixtures
        system = "toy_system"
        models = _write_fixture(fixtures, system)
        run_dir = joinpath(fixtures, system)

        refs = mktempdir()
        ref1 = joinpath(refs, "REFAP_A.pdb")
        ref2 = joinpath(refs, "REFHO_B.pdb")
        cp(joinpath(models, "conf_A1.pdb"), ref1)
        cp(joinpath(models, "conf_C1.pdb"), ref2)

        # Fake DeepAccNet writes only a real results CSV keyed by the flat symlink stem.
        fake_deepaccnet = function (output_dir, sif_path; kwargs...)
            csv_path = joinpath(output_dir, "deepaccnet_results.csv")
            open(csv_path, "w") do io
                println(io, "sample,cb-lddt")
                for group in ("A", "B", "C"), rep = 1:3
                    println(io, "cluster_1_conf_$(group)$(rep),0.9")
                end
            end
            return csv_path
        end

        result = redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                (args...; kwargs...) -> error("prepare must not run"),
                (args...; kwargs...) -> :triaged,
                fake_deepaccnet,
                AlphaConformers.cluster_conformers;
                output_dir = run_dir,
                prepare = false,
                predict = false,
                triage = true,
                query_struct = joinpath(fixtures, "TOY_A.pdb"),
                sifts_uniprot_mapping = DataFrames.DataFrame(
                    PDB = String[],
                    CHAIN = String[],
                ),
                deepaccnet_sif = "image.sif",
                ref1 = ref1,
                ref2 = ref2,
            )
        end

        @test isdir(joinpath(run_dir, "conformer_clustering"))
        @test isfile(
            joinpath(
                run_dir,
                "conformer_clustering",
                "deepaccnet",
                "surviving",
                "surviving.csv",
            ),
        )
        @test result.cluster_result.surviving !== nothing
    end
end

@testitem "alphaconformers runs real DeepAccNet + clustering on a GPU workstation" begin
    import DataFrames

    ENV["GKSwstype"] = "100"  # headless GKS backend so figures render without a display
    # This test drives the real DeepAccNet container and real clustering end to end (the two
    # sub-steps the triage step runs after `triage_outputs`, which is stubbed here so the fixture
    # does not need full triage inputs). It is gated on the local GPU workstation: a GPU,
    # apptainer, the DeepAccNet `.sif` image and the 1AKZ conformer database must all be present,
    # so it auto-skips everywhere else (CI, laptops).
    gpu_ok = try
        success(`nvidia-smi -L`)
    catch
        false
    end
    sif = "/data/MEMBERS/lucas.vitoriano/sif_images/DeepAccNet/deepaccnet.sif"
    data = "/data/MEMBERS/lucas.vitoriano/AlphaConformersDB/1AKZ"
    repo = dirname(@__DIR__)
    ref1 = joinpath(repo, "refs", "1AKZ_A.pdb")
    ref2 = joinpath(repo, "refs", "1SSP_E.pdb")
    resources_ok =
        Sys.which("apptainer") !== nothing &&
        gpu_ok &&
        isfile(sif) &&
        isdir(data) &&
        isfile(ref1) &&
        isfile(ref2)

    # True when a cluster directory holds at least one AlphaFold prediction model.
    function has_models(cluster_dir)
        for inner in readdir(cluster_dir)
            models = joinpath(cluster_dir, inner, "predictions", "sequences", "models")
            isdir(models) &&
                any(f -> endswith(lowercase(f), ".pdb"), readdir(models)) &&
                return true
        end
        return false
    end

    resources_ok && mktempdir() do tmp
        # Copy the populated cluster(s) into the tempdir, preserving the tree, so the real run
        # never mutates the shared database.
        copied = String[]
        for name in readdir(data)
            cluster_dir = joinpath(data, name)
            startswith(name, "cluster_") && isdir(cluster_dir) && has_models(cluster_dir) || continue
            cp(cluster_dir, joinpath(tmp, name))
            push!(copied, name)
        end
        @test !isempty(copied)

        result = AlphaConformers._alphaconformers(
            (args...; kwargs...) -> error("prepare must not run"),
            (args...; kwargs...) -> :triaged,
            AlphaConformers.run_deepaccnet,
            AlphaConformers.cluster_conformers;
            output_dir = tmp,
            prepare = false,
            predict = false,
            triage = true,
            query_struct = ref1,
            sifts_uniprot_mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[]),
            deepaccnet_sif = sif,
            ref1 = ref1,
            ref2 = ref2,
        )

        # DeepAccNet wrote its score table.
        @test isfile(result.deepaccnet_result)
        @test isfile(joinpath(tmp, "deepaccnet", "deepaccnet_results.csv"))

        # Clustering wrote its flat tables under conformer_clustering/.
        results_dir = result.cluster_result.results_dir
        @test isfile(joinpath(results_dir, "aligned_clustering_results.csv"))
        @test isfile(joinpath(results_dir, "aligned_cluster_rmsd.csv"))

        # The DeepAccNet scores auto-fed the score filter, which wrote its surviving table.
        @test isfile(joinpath(results_dir, "deepaccnet", "surviving", "surviving.csv"))
        @test result.cluster_result.surviving !== nothing

        # Both references were assigned to a cluster.
        references = Set(result.cluster_result.reference_clusters.Reference)
        @test "1AKZ_A" in references
        @test "1SSP_E" in references
    end
end
