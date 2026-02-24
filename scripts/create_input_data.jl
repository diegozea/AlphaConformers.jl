#!/stockage/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l ncpus=20
#PBS -l mem=100g
#PBS -l host=node48
#PBS -l walltime=300:00:00
#PBS -j oe
=#

import Pkg
# Pkg.activate(abspath(@__DIR__, ".."))
Pkg.activate("/home/diego.zea/.julia/dev/AlphaConformers/")

import MIToS
import CSV
import DataFrames
import Clustering
import PairwiseListMatrices


# DIR

const WORKING_DIR = "/alpha/runs/diego.zea/from_apo"

cd(WORKING_DIR)

# LOAD DATA

const DATA = DataFrames.DataFrame(CSV.File("/home/diego.zea/.julia/dev/AlphaConformers/data/Supplementary_Table_1_91_apo_holo_pairs.csv", delim=';'))

download("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz", 
    "pdb_chain_uniprot.csv.gz")
const sifts_up_data = CSV.File("pdb_chain_uniprot.csv.gz", comment="#") |> DataFrames.DataFrame


# FUNCTIONS

function run_foldseek_search(pdb_file)
    foldseek = "/stockage/EQUIPES/AMIG/PROGRAMMES/foldseek/bin/foldseek"
    # db = "/stockage/EQUIPES/AMIG/MEMBERS/carla.martins/foldseek/db_pdb/pdb"
    db = "/alpha/database/pdb"
    tmp_folder = mktempdir()
    out_file = "$(pdb_file)_results.m8"
    run(`$foldseek easy-search $pdb_file $db $out_file $tmp_folder`)
    out_file
end

"""
    parse_tm_align_output(tm_align_output::String)

Parse the output of a TM-align alignment using `-dir` and return a named tuple containing
the name of Chain_1, the name of Chain_2, RMSD, Seq_ID, and the TM-scores 
(normalized by the length of Chain_1 and Chain_2) for each alignment.

# Arguments
- `tm_align_output::String`: The text output from TM-align.

# Returns
- `Vector{NamedTuple}`: A vector of named tuples, where each tuple corresponds to 
  an alignment and contains the following fields: `:chain_a`, `:chain_b`, `:RMSD`, 
  `:seq_id`, `:TM_score_a`, and `:TM_score_b`.
"""
function parse_tm_align_output(tm_align_output::String)
    alignments = split(tm_align_output, r" \*+") 
    # TODO: improve this by iterating the lines and looking for Chain_1 as an alignment start
    parsed_alignments = NamedTuple{(:chain_a, :chain_b, :RMSD, :seq_id, :TM_score_a, :TM_score_b), Tuple{String, String, Vararg{Float64, 4}}}[]
    
    for alignment in alignments
        if alignment != "" && occursin("Chain_1", alignment)
            chain_a = match(r"Chain_1: .*/(.*\.pdb_?[a-zA-Z0-9]?)", alignment).captures[1]
            chain_b = match(r"Chain_2: .*/(.*\.pdb_?[a-zA-Z0-9]?)", alignment).captures[1]
            RMSD = parse(Float64, match(r"RMSD=\s+(.*),", alignment).captures[1])
            seq_id = parse(Float64, match(r"Seq_ID=n_identical/n_aligned=\s+(.*)", alignment).captures[1])
            TM_score_a = parse(Float64, match(r"TM-score=\s+(.*) \(if normalized by length of Chain_1", alignment).captures[1])
            TM_score_b = parse(Float64, match(r"TM-score=\s+(.*) \(if normalized by length of Chain_2", alignment).captures[1])
            push!(parsed_alignments, (chain_a=chain_a, chain_b=chain_b, RMSD=RMSD, seq_id=seq_id, TM_score_a=TM_score_a, TM_score_b=TM_score_b))
        end
    end
    
    return parsed_alignments
end

# create a function to look up the UniProt accession for a given PDB code and chain (if available)
function get_uniprot_acc(pdb_file::String, data::DataFrames.DataFrame=sifts_up_data)
    pdb_code = lowercase(match(r"([a-zA-Z0-9]{4})\.pdb", pdb_file).captures[1])
    chain_match = match(r".*\.pdb_([a-zA-Z0-9])", pdb_file)
    if chain_match !== nothing
        chain_code = chain_match.captures[1]
        ups = data[(data.PDB .== pdb_code) .& (data.CHAIN .== chain_code), :SP_PRIMARY]
    else
        ups = data[data.PDB .== pdb_code, :SP_PRIMARY]
    end
    unique(ups)
end

# create a function to look up the structures associated with a given UniProt accession
function get_pdb_codes(uni_acc, data::DataFrames.DataFrame=sifts_up_data)
    pdb_codes = data[data.SP_PRIMARY .== uni_acc, [:PDB, :CHAIN]]
    uppercase.(unique(pdb_codes))
end


# delete the query if it is in the target column
function delete_query_from_target!(search_results, pdb_code, chain_code)
    up = first(get_uniprot_acc("$(pdb_code).pdb_$(chain_code)"))
    query_structures = get_pdb_codes(up)
    for row in eachrow(query_structures)
        pdb = "$(row[:PDB]).pdb"
        pdb_chain = "$(row[:PDB]).pdb_$(row[:CHAIN])"
        filter!(row -> row.target != pdb && row.target != pdb_chain, search_results)
    end
    search_results
end

# look for known conformations of the proteins showing similar structures to the query protein
function list_known_conformations(search_results)
    uniprots = Set{String}()
    new_targets = Set{String}()
    for row in eachrow(search_results)
        pdb_code = match(r"([a-zA-Z0-9]{4})\.pdb", row.target).captures[1]
        match_chain = match(r".*\.pdb_([a-zA-Z0-9])", row.target)
        if match_chain !== nothing
            chain_code = match_chain.captures[1]
            target = "$(pdb_code).pdb_$(chain_code)"
        else
            target = "$(pdb_code).pdb"
        end
        ups = get_uniprot_acc(target)
        for up in ups
            push!(uniprots, up)
        end
    end
    for up in uniprots
        pdbs = get_pdb_codes(up)
        for row in eachrow(pdbs)
            pdb = row[:PDB]
            chain = row[:CHAIN]
            target = "$(pdb).pdb"
            target_with_chain = "$(pdb).pdb_$(chain)"
            if !(target in search_results.target) && !(target_with_chain in search_results.target)
                push!(new_targets, target_with_chain)
            end
        end
    end
    new_targets
end


function create_pdb_folder(targets::Set{String})
    # PATH to the PDB files
    pdb_folder = "/alpha/database/pdb/pdb_files"

    # Create a local folder for PDB files 
    # NOTE : It changes from first_targets to targets here
    if isdir("targets")
        rm("targets"; recursive=true)
    end
    mkdir("targets")
    local_pdb_folder = abspath("targets")
    
    for target in targets
        # Split target into pdb_code and chain_code
        fields = split(target, '_')
        pdb_code = String(fields[1])
        if length(fields) == 2
            chain_code = String(fields[2])
        else
            chain_code = MIToS.PDB.All
        end
        # Define original and local file paths
        original_pdb_path = joinpath(pdb_folder, pdb_code)
        if !isfile(original_pdb_path)
            continue
        end
        local_pdb_path = joinpath(local_pdb_folder, target)
        # Read only specific chain and write to local file
        # occupancyfilter is needed to avoid the duplicated residue warnings with TMalign
        chain = MIToS.PDB.read(original_pdb_path, MIToS.PDB.PDBFile, 
            onlyheavy=true, occupancyfilter=true,
            chain=chain_code)
        MIToS.PDB.write(local_pdb_path, chain, MIToS.PDB.PDBFile)
    end

    return local_pdb_folder
end


function run_tmalign_and_cluster(local_pdb_folder::String, filename_prefix::String="tmalign")
    # create a temporary file for the chain list
    chain_list_path = tempname()
    open(chain_list_path, "w") do io
        for file in readdir(local_pdb_folder)
            if occursin(".pdb", file) # only select PDB files
                println(io, file)
            end
        end
    end

    # run TMalign with the -dir option:
    dir_path = abspath(local_pdb_folder)
    # add / at the end of the path if it is not there
    if dir_path[end] != '/'
        dir_path *= '/'
    end
    
    stdout_path = filename_prefix * ".out"
    stderr_path = filename_prefix * ".err"
    
    run(pipeline(`/stockage/EQUIPES/AMIG/MEMBERS/diego.zea/bin/TMalign -dir $dir_path $chain_list_path -fast`, 
        stdout=stdout_path, stderr=stderr_path))
    tmalign_out = read(stdout_path, String)
    
    tmalign_data = DataFrames.DataFrame(parse_tm_align_output(tmalign_out))

    rmsd_mat = PairwiseListMatrices.from_table(tmalign_data[:,1:3], false)

    rmsd_clust = Clustering.hclust(rmsd_mat, linkage=:complete, branchorder=:optimal)

    rmsd_clusters = Clustering.cutree(rmsd_clust, h=1.0)

    rmsd_chains = names(rmsd_mat)[1]
    
    # Return named tuple
    return (tmalign_data = tmalign_data, rmsd_mat = rmsd_mat, rmsd_chains = rmsd_chains, rmsd_clusters = rmsd_clusters)
end


function main(row)
    apo_id = row[:apo_id]
    m = match(r"([0-9A-Z]{4}).*_([A-Z])", apo_id)
    if m !== nothing
        pdb_code = String(m.captures[1])
        chain_code = String(m.captures[2])
    else
        return nothing
    end

    folder = "$(pdb_code)_$(chain_code)"
    if isdir(folder)
        rm(folder; recursive=true)
    end
    mkdir(folder)
    cd(folder)
    try
        # download pdb file
        pdb_file = MIToS.PDB.downloadpdb(pdb_code, format=MIToS.PDB.PDBFile)
        chain = read("$(pdb_code).pdb.gz", MIToS.PDB.PDBFile, chain=chain_code)
        pdb_file = abspath("$(pdb_code)_$(chain_code).pdb")
        write("$(pdb_code)_$(chain_code).pdb", chain, MIToS.PDB.PDBFile)

        # run foldseek locally
        search_results_file = run_foldseek_search(pdb_file)

        # parse foldseek results
        column_names = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits"]
        search_results = DataFrames.DataFrame(CSV.File(search_results_file, delim='\t', header=column_names))

        # delete query from target
        delete_query_from_target!(search_results, pdb_code, chain_code)

        # download new targets
        new_targets = list_known_conformations(search_results)

        # create local folder for new targets
        targets = union(Set{String}(search_results.target), new_targets)
        local_pdb_folder = create_pdb_folder(targets)

        # run TMalign and cluster
        tmalign_results = run_tmalign_and_cluster(local_pdb_folder)

        # get query sequence
        query_res = read(pdb_file, MIToS.PDB.PDBFile, group="ATOM", onlyheavy=true, occupancyfilter=true)
        query_seq = join(MIToS.MSA.three2residue(residue.id.name) for residue in query_res)

        # create clusters folder
        if isdir("clusters")
            rm("clusters"; recursive=true)
        end
        mkdir("clusters")

        # create the input data to run alphafold for each cluster
        for cluster in unique(tmalign_results.rmsd_clusters)
            targets = tmalign_results.rmsd_chains[tmalign_results.rmsd_clusters .== cluster]
            cluster_folder = "clusters/cluster_$(cluster)"
            mkdir(cluster_folder)
            open(joinpath(cluster_folder, "sequences.fasta"), "w") do file 
                println(file, ">query")
                println(file, query_seq)
                for target in targets
                    res = read(joinpath(local_pdb_folder, target), MIToS.PDB.PDBFile, group="ATOM", onlyheavy=true, occupancyfilter=true)
                    println(file, ">$target")
                    for residue in res
                        print(file, MIToS.MSA.three2residue(residue.id.name))
                    end
                    print(file, "\n")
                end
            end
            run(`clustalo -i $(joinpath(cluster_folder, "sequences.fasta")) -o $(joinpath(cluster_folder, "sequences.aln")) --force --output-order=input-order`)
            msa = read(joinpath(cluster_folder, "sequences.aln"), MIToS.MSA.FASTA)
            MIToS.MSA.adjustreference!(msa)
            complete_msa = msa[vec(MIToS.MSA.coverage(msa) .> 0.5), :]
            MIToS.MSA.write(joinpath(cluster_folder, "sequences.a3m"), complete_msa, MIToS.MSA.FASTA)
            mkdir(joinpath(cluster_folder, "templates"))
            pdbs = unique(String["$(split(pdb, '.')[1]).pdb" for pdb in targets])
            for pdb in pdbs
                res = read(joinpath("/alpha/database/pdb/pdb_files", pdb), MIToS.PDB.PDBFile, model="1", onlyheavy=true, occupancyfilter=true)
                MIToS.PDB.write(joinpath(cluster_folder, "templates", lowercase(pdb)), res, MIToS.PDB.PDBFile)
            end
        end
    catch err
        @warn err
    finally
        cd("..")
    end
end

# RUN

for row in DataFrames.eachrow(DATA)
    main(row)
end