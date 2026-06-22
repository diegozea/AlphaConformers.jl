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

    mktempdir() do dir
        pdb = joinpath(dir, "numbered.pdb")
        write_numbered_pdb(pdb)
        @test AlphaConformers.get_res_beg_end(pdb) == (1, 10)
    end
end
