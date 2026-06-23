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
        return fake_preparer, fake_predictor, fake_triager
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

        for combination in combinations
            calls = Dict{Symbol,Any}()
            fake_preparer, fake_predictor, fake_triager = fake_steps(calls)

            result = redirect_stdout(devnull) do
                if combination.predict
                    AlphaConformers._alphaconformers(
                        fake_preparer,
                        fake_triager,
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
                        predictor_keyword = 42,
                    )
                else
                    AlphaConformers._alphaconformers(
                        fake_preparer,
                        fake_triager;
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
                    )
                end
            end

            @test result.output_dir == abspath(output_dir)
            @test result.prepared == combination.prepare
            @test result.predicted == combination.predict
            @test (result.triage_result !== nothing) == combination.triage
            @test haskey(calls, :prepare) == combination.prepare
            @test haskey(calls, :predict) == combination.predict
            @test haskey(calls, :triage) == combination.triage

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

    mktempdir() do dir
        query_struct = joinpath(dir, "1ABC_A.pdb")
        pdb_folder = joinpath(dir, "pdb")
        output_dir = joinpath(dir, "run")
        mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[])

        result = redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                fake_preparer,
                fake_triager,
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

    mktempdir() do dir
        query_struct = joinpath(dir, "1ABC_A.pdb")
        pdb_folder = joinpath(dir, "pdb")
        output_dir = joinpath(dir, "run")
        mapping = DataFrames.DataFrame(PDB = String[], CHAIN = String[])

        result = redirect_stdout(devnull) do
            AlphaConformers._alphaconformers(
                fake_preparer,
                fake_triager,
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
            )
        end
        @test prepare_called[] == false
        @test occursin("foldseek_results_folder", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                fail_if_prepare_runs,
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
            )
        end
        @test prepare_called[] == false
        @test occursin("foldseek_results_folder", msg)

        msg = argument_error_message() do
            AlphaConformers._alphaconformers(
                fail_if_prepare_runs,
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
            )
        end
        @test prepare_called[] == false
        @test occursin("foldseek_results_folder", msg)
    end
end
