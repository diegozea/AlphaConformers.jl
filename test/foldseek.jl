using TestItems

@testitem "Foldseek execution against generated database" begin
    import DataFrames
    import Foldseek_jll
    import Printf

    function write_test_pdb(path)
        residues = [
            "ALA",
            "GLY",
            "SER",
            "THR",
            "LEU",
            "VAL",
            "ASP",
            "LYS",
            "GLU",
            "ILE",
            "ASN",
            "PHE",
            "TYR",
            "ARG",
            "GLN",
        ]
        serial = 1
        open(path, "w") do io
            for (resnum, resname) in enumerate(residues)
                x = 3.8 * (resnum - 1)
                atoms = (
                    ("N", x - 1.2, 0.0, 0.0, "N"),
                    ("CA", x, 0.0, 0.0, "C"),
                    ("C", x + 1.2, 1.0, 0.0, "C"),
                    ("O", x + 1.8, 2.0, 0.0, "O"),
                )
                for (atom, ax, ay, az, element) in atoms
                    Printf.@printf(
                        io,
                        "ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n",
                        serial,
                        atom,
                        resname,
                        resnum,
                        ax,
                        ay,
                        az,
                        element,
                    )
                    serial += 1
                end
            end
            println(io, "TER")
            println(io, "END")
        end
    end

    mktempdir() do dir
        query_pdb = joinpath(dir, "query.pdb")
        target_pdb = joinpath(dir, "target.pdb")
        db_path = joinpath(dir, "test_db")
        db_path2 = joinpath(dir, "test_db2")
        write_test_pdb(query_pdb)
        cp(query_pdb, target_pdb)
        run(`$(Foldseek_jll.foldseek()) createdb $target_pdb $db_path`)
        run(`$(Foldseek_jll.foldseek()) createdb $target_pdb $db_path2`)

        search_output = foldseek_search(
            query_pdb;
            db_path = db_path,
            format_output = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qtmscore,ttmscore,alntmscore,rmsd,prob",
        )
        @test search_output == joinpath(dir, "query_results.m8")
        @test isfile(search_output)
        @test filesize(search_output) > 0
        @test DataFrames.nrow(read_foldseek_search_results(search_output)) >= 1

        run_output = run_foldseek(query_pdb, 1, db_path; out_folder = dir)
        @test length(run_output) == 1
        @test isfile(run_output[1].table_file)
        @test filesize(run_output[1].table_file) > 0
        @test isfile(run_output[1].msa_file)
        @test isdir(run_output[1].aligned_structures_folder)
        @test !isempty(readdir(run_output[1].aligned_structures_folder))

        nofilter_dir = joinpath(dir, "nofilter")
        mkpath(nofilter_dir)
        nofilter_output =
            run_foldseek(query_pdb, 1, db_path; out_folder = nofilter_dir, filtrage = false)
        @test length(nofilter_output) == 1
        @test isfile(nofilter_output[1].table_file)

        vector_dir = joinpath(dir, "vector")
        mkpath(vector_dir)
        vector_output = run_foldseek(query_pdb, 1, [db_path]; out_folder = vector_dir)
        @test length(vector_output) == 1
        @test isfile(vector_output[1].table_file)

        cd(dir) do
            relative_output = run_foldseek(
                "query.pdb",
                1,
                "test_db";
                out_folder = "missing_parent/relative",
            )
            @test length(relative_output) == 1
            @test isfile(relative_output[1].table_file)
            @test isdir(relative_output[1].aligned_structures_folder)
            @test realpath(dirname(dirname(relative_output[1].table_file))) ==
                  realpath(joinpath(dir, "missing_parent", "relative"))

            relative_vector_output = run_foldseek(
                "query.pdb",
                1,
                ["test_db"];
                out_folder = "missing_parent/vector_relative",
            )
            @test length(relative_vector_output) == 1
            @test isfile(relative_vector_output[1].table_file)

            relative_comma_output = run_foldseek(
                "query.pdb",
                1,
                "test_db, test_db2";
                out_folder = "missing_parent/comma_relative",
            )
            @test length(relative_comma_output) == 2
            @test all(output -> isfile(output.table_file), relative_comma_output)
        end
    end
end

@testitem "Foldseek table parsing and merging" begin
    import DataFrames

    function write_m8(path, rows)
        open(path, "w") do io
            for row in rows
                println(io, join(row, '\t'))
            end
        end
    end

    row_a = [
        "query.pdb",
        "target_a.pdb_A",
        0.95,
        100,
        1,
        0,
        1,
        100,
        5,
        104,
        1e-20,
        250,
        0.8,
        0.7,
        0.75,
        1.2,
        0.99,
    ]
    row_b = [
        "query.pdb",
        "target_b.pdb_A",
        0.75,
        80,
        20,
        1,
        3,
        82,
        10,
        89,
        1e-4,
        100,
        0.5,
        0.4,
        0.45,
        3.0,
        0.8,
    ]

    mktempdir() do dir
        table_a = joinpath(dir, "a.m8")
        table_b = joinpath(dir, "b.m8")
        write_m8(table_a, [row_a, row_b])
        write_m8(table_b, [row_a])

        parsed = read_foldseek_search_results(table_a)
        @test DataFrames.nrow(parsed) == 2
        @test names(parsed) == [
            "query",
            "target",
            "fident",
            "alnlen",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "evalue",
            "bits",
            "qtmscore",
            "ttmscore",
            "alntmscore",
            "rmsd",
            "prob",
        ]
        @test parsed[1, :target] == "target_a.pdb_A"
        @test parsed[1, :rmsd] == 1.2

        merged = merge_tables([table_a, table_b])
        @test DataFrames.nrow(merged) == 2
        @test DataFrames.ncol(merged) == 18
        @test Set(merged.target) == Set(["target_a.pdb_A", "target_b.pdb_A"])
        @test only(unique(merged.file[merged.target .== "target_a.pdb_A"])) ==
              abspath(table_a)
    end
end

@testitem "Foldseek argument validation and alignment positions" begin
    mktempdir() do dir
        pdb_file = joinpath(dir, "query.pdb")
        db_dir = joinpath(dir, "db")
        touch(pdb_file)
        mkdir(db_dir)

        @test_throws ErrorException foldseek_search(pdb_file)
        @test_throws ErrorException foldseek_search(pdb_file; db_path = db_dir)
        @test_throws ErrorException run_foldseek(pdb_file, 1, db_dir)
    end

    aln = [('A', 'A'), ('C', '-'), ('G', 'G'), ('-', 'T'), ('T', 'T')]
    @test get_aligned_positions(aln) == [(1, 1), (3, 2), (4, 4)]
end
