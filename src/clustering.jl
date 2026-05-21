"""
    coverage_and_identity(aln) -> (coverage, identity)

Compute the coverage and sequence identity of an pairwise alignment,
from the perspective of sequence `a` (the query).  Gaps in `a` are ignored entirely, because we want to measure 
how much of the *query* is covered by the subject `b`, not the reverse.
Input : 
- `aln`: an iterable of `(res_a, res_b)` residue pairs representing a pairwise alignment.
Output : 
A tuple `(coverage, identity)` where:
- `coverage`  : fraction of query residues (non-gap in `a`) that are matched to a non-gap residue in `b`.
- `identity`  : fraction of **matched positions** that are identical.
Interpret result 
- Low coverage + high identity → `b` aligns well but only on a small region of `a`.
- High coverage + low identity → `b` spans most of `a` but with many mismatches.
"""
function coverage_and_identity(aln)
    total_a = 0
    matched_a = 0
    identities = 0
    for (res_a, res_b) in aln
        if res_a != '-'
            total_a += 1
            if res_b != '-'
                matched_a += 1
                if res_a == res_b
                    identities += 1
                end
            end
        end
    end
    (matched_a / total_a, identities / matched_a)
end
"""
    structural_alignment(conformation_a, conformation_b; aln_type, aln_model, limit_residues_range)
    -> (aligned_a, aligned_b, matches, rmsd, coverage, identity) or nothing

Align two protein structures by sequence, then superimpose the matched residues.
Input : 
- `conformation_a` : query structure (vector of residues, e.g. from MIToS.PDB).
- `conformation_b` : target structure to align onto `a`.
Output : 
A named tuple `(aligned_a, aligned_b, matches, rmsd, coverage, identity)` where:
- `aligned_a`, `aligned_b` : superimposed residue arrays.
- `matches`                : vector of `(index_in_a, index_in_b)` aligned residue pairs.
- `rmsd`                   : RMSD after superposition (Å).
- `coverage`, `identity`   : as defined in `coverage_and_identity` (query-relative).
Returns `nothing` if either structure has no Cα atoms, sequences cannot be extracted,
no residues fall within `limit_residues_range`, or an unexpected error occurs.
"""
function structural_alignment(conformation_a, conformation_b;
    aln_type=BioAlignments.OverlapAlignment(),
    aln_model=BioAlignments.AffineGapScoreModel(match=6, mismatch=-4, gap_open=-2,
        gap_extend=-1),
    limit_residues_range::Union{UnitRange{Int}, Tuple{Int,Int}, Nothing}=nothing)
    
    # only keep residues with 'CA' atom
    clean_a = filter(res -> !isempty(MIToS.PDB.findatoms(res, "CA")), conformation_a)
    clean_b = filter(res -> !isempty(MIToS.PDB.findatoms(res, "CA")), conformation_b)
    len_a = length(clean_a)
    len_b = length(clean_b)
    if len_a == 0 || len_b == 0
        @warn "One of the structures has no 'CA' atoms: $len_a, $len_b"
        return nothing
    end
    
    try
        # get the sequences
        seqs_a = MIToS.PDB.modelled_sequences(clean_a)
        seqs_b = MIToS.PDB.modelled_sequences(clean_b)

        if isempty(seqs_a) || isempty(seqs_b)
            @warn "No modelled sequences could be extracted (empty sequences)"
            return nothing
        end

        seq_a = first(values(seqs_a))
        seq_b = first(values(seqs_b))

        # align the sequences
        aln = BioAlignments.alignment(BioAlignments.pairalign(aln_type, seq_a, seq_b, aln_model))
        # get the stats
        coverage, identity = coverage_and_identity(aln)
        # get the aligned residues
        matches = get_aligned_positions(aln)

        # filter matches to the requested residue range (positions in conformation_a)
        if limit_residues_range !== nothing
            start_res, stop_res = if limit_residues_range isa UnitRange{Int}
                first(limit_residues_range), last(limit_residues_range)
            else
                limit_residues_range  # Tuple{Int,Int}
                
            end

            if start_res < 1 || stop_res > len_a || start_res > stop_res
                @warn "Invalid residue range ($start_res, $stop_res) for structure of length $len_a — ignoring range filter"
            else
                # matches is a vector of (index_in_a, index_in_b) pairs
                matches = filter(m -> first(m) >= start_res && first(m) <= stop_res, matches)
                
                if isempty(matches)
                    @warn "No aligned residues found in range ($start_res, $stop_res)"
                    return nothing
                end
            end
        end

        # structural superposition of the aligned residues
        aligned_a, aligned_b, rmsd = MIToS.PDB.superimpose(clean_a, clean_b, matches)
        return (aligned_a, aligned_b, matches, rmsd, coverage, identity)
    catch err
        @error "Error in the structural alignment: $err"
        return nothing
    end
end

function within_cluster(item1,item2,threshold)
    _,_,_,rmsd,coverage,_=AlphaConformers.structural_alignment(item1, item2)
    return rmsd <= threshold
end

function get_target2sequence(expanded_table, msa)
    seqnames=collect(MIToS.MSA.sequencename_iterator(msa))
    clean_seqnames = collect(split(x, '\t')[1] for x in seqnames)

    target2sequence = Dict{String,String}()
    test=[]
    for row in eachrow(expanded_table)
        seqname = ismissing(row.evalue) ? row.query : row.target
        seqname = String(seqname)
        push!(test,seqname)
        if seqname in clean_seqnames
            target2sequence[row.target] = seqname
        else 
            @warn "The sequence $seqname is not in the MSA."
        end
    end
    target2sequence
    
end

function get_cluster2targets(targets, clusters)
    cluster2targets = OrderedCollections.OrderedDict{Int,Vector{String}}()
    cluster_numbers = unique(clusters)
    for cluster in cluster_numbers
        cluster2targets[cluster] = String[]
        for (name, cl) in targets
            if cl == cluster
                if name in cluster2targets[cluster]
                    continue
                end
                push!(cluster2targets[cluster], name)
            end
        end
    end

    cluster2targets, length(cluster2targets)

end

function get_cluster2targets_concatenate(targets, clusters)
    cluster2targets = OrderedCollections.OrderedDict{Int,Vector{String}}()
    # Regroupe les clusters par paquets de 10
    grouped_clusters = Dict{Int, Vector{String}}()
    for (name, cl) in targets
        group = Int(ceil(cl / 10))
        if !haskey(grouped_clusters, group)
            grouped_clusters[group] = String[]
        end
        # Ajoute sans doublon
        if !(name in grouped_clusters[group])
            push!(grouped_clusters[group], name)
        end
    end
    # Convertit en OrderedDict pour compatibilité
    for (group, names) in sort(collect(grouped_clusters))
        cluster2targets[group] = names
    end
    @show length(cluster2targets)
    @show first(cluster2targets, 5)
    cluster2targets, length(cluster2targets)
end

function get_cluster2seqnames(cluster2targets, target2sequence)
    cluster2seqnames = OrderedCollections.OrderedDict{Int,Vector{String}}()
    for (cluster, targets) in cluster2targets
        cluster2seqnames[cluster] = unique(target2sequence[target] for target in targets
                                           if haskey(target2sequence, target))
    end
    cluster2seqnames
end

function _check_names_in_msa(names, msa)
    found_names = []
    seqnames= collect(MIToS.MSA.sequencename_iterator(msa))
    clean_seqnames = collect(split(x, '\t')[1] for x in seqnames)
    for seqname in clean_seqnames
        if seqname in names
            push!(found_names, seqname)
        end
    end
    missing_names = setdiff(names, found_names)
    if !isempty(missing_names)
        @info "The first 5 sequences in the MSA are: $(first(MIToS.MSA.sequencename_iterator(msa), 5))"
        @warn "The following sequences are not in the MSA: $(collect(missing_names))"
    end
    isempty(missing_names)
end

function _rename_msa!(msa)
    new_names = MIToS.MSA.sequencenames(msa)
    MIToS.MSA.rename_sequences!(msa, new_names)
end

function get_cluster2msa(msa, cluster2seqnames)
    query_name = first(MIToS.MSA.sequencename_iterator(msa))
    cluster2msa = OrderedCollections.OrderedDict{Int,MIToS.MSA.AnnotatedMultipleSequenceAlignment}()
    for (cluster, seqnames) in cluster2seqnames
        seqnames = copy(seqnames)
        pushfirst!(seqnames, query_name)
        _rename_msa!(msa)
        if ! _check_names_in_msa(seqnames, msa)
            @warn "The cluster $cluster will be skipped."
            continue
        end
        try
            names=MIToS.MSA.sequencenames(msa)
            all_name= String[]
            for nom in seqnames
                name_line = findfirst(x -> startswith(x, nom), names)
                if name_line !== nothing
                    push!(all_name, names[name_line])
                end
            end
            if isempty(all_name)
                @warn "No MSA for cluster $cluster"
                continue
            end
            all_name = unique(all_name)
            cluster2msa[cluster] = msa[all_name, :]
        catch err
            @error "Error ($err) getting the MSA for cluster $cluster"
            @info "seqnames: $seqnames"
            @info "msa names: $(MIToS.MSA.sequencenames(msa))"
            rethrow(err)
        end
    end
    cluster2msa
end
"""
We can only have 4 templates for each cluster to run AlphaFold with ColabFold
To select those templates we divide the cluster in 4 equal part and select one template for each part to have a good representation of the cluster
If less than 8 templates are available we take the four first one, if less than 4 we take all the templates
"""
#Select 4 template for each clusters
function get_cluster2structures(structures, cluster2targets)
    # dictionnaire final
    cluster2structures = OrderedCollections.OrderedDict{Int, Dict{String, Vector{MIToS.PDB.PDBResidue}}}()

    # pour tracker les targets déjà assignés
    used_targets = Set{String}()

    for (cluster, targets) in cluster2targets
        # on filtre les targets déjà utilisés
        available_targets = filter(t -> !(t in used_targets), targets)
        
        # on choisit les sous-targets
        subtargets = String[]
        n = length(available_targets)

        if n == 0
            @warn "Cluster $cluster has no available targets left after filtering"
            continue
        elseif n < 8
            subtargets = available_targets[1:min(4, n)]
        else
            step = max(1, div(n, 4))
            for i in 1:4
                idx = i * step
                if idx <= n
                    push!(subtargets, available_targets[idx])
                end
            end
        end

        # marquer ces targets comme utilisés
        foreach(t -> push!(used_targets, t), subtargets)

        # construire cluster2structures
        cluster2structures[cluster] = Dict(
            t => structures[t] for t in subtargets
        )
    end

    cluster2structures
end

"""
    create_template_clusters_hobohm(expanded_table, msa, structures, cutoff)
    -> (nb_cluster, cl2msa, cl2pdb)

Cluster structural templates using the Hobohm I algorithm and build per-cluster
MSAs and structure dictionaries.
Input : 
- `expanded_table` : DataFrame of Foldseek hits with at least a `target` column.
- `msa`            : full MSA from which per-cluster sub-MSAs will be extracted.
- `structures`     : ordered dict mapping target names to their residue arrays.
- `cutoff`         : similarity cutoff for Hobohm I clustering.
Output : 
- `nb_cluster` : total number of clusters found.
- `cl2msa`     : dict mapping cluster index → sub-MSA (query + cluster members).
- `cl2pdb`     : dict mapping cluster index → dict of target name → residue array.
"""
function create_template_clusters_hobohm(
    expanded_table::DataFrames.DataFrame,
    msa::MIToS.MSA.AnnotatedMultipleSequenceAlignment,
    structures::OrderedCollections.OrderedDict{String,Vector{MIToS.PDB.PDBResidue}},
    cutoff::Float64)
    
    target2sequence = AlphaConformers.get_target2sequence(expanded_table, msa)
    targets = Set{String}(expanded_table.target)
    
    #Hobohm I clustering
    println("Clustering structures with Hobohm I algorithm...")
    names = collect(keys(structures))
    struct_list = collect(values(structures))
    output=MIToS.MSA.hobohmI(within_cluster, struct_list, cutoff,threads=true)
    assignments = output.clusters
    cluster_labels = Dict{String, Int}()

    for (name, cluster_id) in zip(names, assignments)
        cluster_labels[name] = cluster_id
    end

    clusters = [cluster_labels[t] for t in targets]
    println("Number of clusters found: $(length(unique(clusters)))")

    cluster2targets, nb_cluster = AlphaConformers.get_cluster2targets(cluster_labels, clusters)
    
    cl2seq = AlphaConformers.get_cluster2seqnames(cluster2targets, target2sequence)
    println("Create the related MSAs for each cluster...")
    cl2msa = AlphaConformers.get_cluster2msa(msa, cl2seq)
    println("Getting the templates for each cluster...")
    cl2pdb = AlphaConformers.get_cluster2structures(structures, cluster2targets)
    (nb_cluster, cl2msa, cl2pdb)
end

"""
    create_folder_structure_hobohm(clusters, cl2msa, cl2pdb; out_folder)
    -> out_folder

Write per-cluster MSA and template files to disk, downloading full structure files
from RCSB or AlphaFold DB when available.
Input : 
- `clusters` : total number of clusters (iterates from 1 to `clusters`).
- `cl2msa`   : dict mapping cluster index → MSA (output of `create_template_clusters_hobohm`).
- `cl2pdb`   : dict mapping cluster index → dict of target name → residue array.
- `out_folder` : output directory (default: a temporary folder).
Output : 
out_folder/
├── calpha_template.csv     ← list of templates written from Cα-only coordinates
│                             (i.e. for which the full file download failed)
└── cluster_N/
    ├── sequences.a3m
    └── templates_complete/
        └── <template>.cif or .pdb
"""
function create_folder_structure_hobohm(clusters,
    cl2msa::OrderedCollections.OrderedDict{Int,MIToS.MSA.AnnotatedMultipleSequenceAlignment},
    cl2pdb::OrderedCollections.OrderedDict{Int,Dict{String,Vector{MIToS.PDB.PDBResidue}}};
    out_folder::String=mktempdir())
    #unique_cluster=unique(clusters)
    calpha_template = String[]

    for clust in 1:clusters
        # MSA
        
        cluster_folder = mkdir(joinpath(out_folder, "cluster_$(clust)"))
        msa_file = joinpath(cluster_folder, "sequences.a3m")
        MIToS.MSA.write_file(msa_file, cl2msa[clust], MIToS.MSA.FASTA)

        # Structures
        cluster_template_folder = mkdir(joinpath(cluster_folder, "templates_complete"))
        for (target, structure) in cl2pdb[clust]
            #Get the right extension 
            base_name=first(splitext(basename(target)))
            if startswith(basename(base_name),"AF") 
                name_cif_file=basename(base_name)*".pdb"
                try 
                    MIToS.PDB.download_alphafold_structure(name_cif_file, joinpath(cluster_template_folder, name_cif_file))
                catch e
                    @warn "Problem to download for AFDB the complete template file: "*base_name
                    out_cif = joinpath(cluster_template_folder, name_cif_file)
                    MIToS.PDB.write_file(out_cif, structure, MIToS.PDB.PDBFile)
                    push!(calpha_template, out_cif)

                    continue
                end
                
            else 
                pdbcode=first(splitext(basename(base_name)))
                pdbcode = split(basename(pdbcode), '_')[1]
                name_cif_file=pdbcode*".cif"
                try 
                    MIToS.Utils.download_file("https://files.rcsb.org/download/"*name_cif_file, joinpath(cluster_template_folder,name_cif_file))

                catch e
                    @warn "Problem to download from PDB the template file: "*base_name
                    out_cif = joinpath(cluster_template_folder, name_cif_file)
                    MIToS.PDB.write_file(out_cif, structure, MIToS.PDB.MMCIFFile)
                    push!(calpha_template, out_cif)
                    continue
                end
                
            end
        end
        
    end
    CSV.write(joinpath(out_folder, "calpha_template.csv"),
              (; file = calpha_template))
    out_folder
end

function write_new_msa(new_msa_path,query_id,query_seq,cluster_ids,cluster_seqs)
    open(new_msa_path, "w") do io
        println(io, ">$query_id")
        println(io, query_seq)
        for (id, seq) in zip(cluster_ids, cluster_seqs)
            println(io, ">$id")
            println(io, seq)
        end
    end
end




function download_template(cluster_folder,templates)
    cluster_template_folder = mkdir(joinpath(cluster_folder, "templates_complete"))
    for target in templates
        #Get the right extension 
        base_name = first(splitext(basename(target)))
        if startswith(basename(base_name),"AF") 
            new_name= replace(basename(base_name), "MODEL_V6" => "model_v6")
            name_cif_file=basename(new_name)*".pdb"
            try 
                MIToS.PDB.download_alphafold_structure(name_cif_file, joinpath(cluster_template_folder, name_cif_file))
            catch e
                @warn "Problem to download for AFDB the complete template file: "*name_cif_file
                continue
            end
        else 
            pdbcode= first(splitext(basename(base_name)))
            pdbcode = split(basename(pdbcode), '_')[1]
            @show pdbcode
            name_cif_file=pdbcode*".cif"
            try 
                MIToS.Utils.download_file("https://files.rcsb.org/download/"*name_cif_file, joinpath(cluster_template_folder,name_cif_file))

            catch e
                @warn "Problem to download from PDB the template file: "*name_cif_file
                continue
            end
            
        end
    end
    
end

#Clean the file name in the MSA and prep the template file for AlphaFold (.cif + good name)
"""
    clean_msa_template_names(clusters, out_folder; cluster_name)

Normalize sequence and template names across MSA and template files for all clusters,
ensuring compatibility with downstream structure prediction tools (AlphaFold3, Boltz2).
Input : 
- `clusters`   : number of clusters to process (iterates from 1 to `clusters`).
- `out_folder` : path to the folder containing `cluster_1/`, `cluster_2/`, ... subfolders.
- `cluster_name` : name of the output subfolder for cleaned templates
                   (default: `"templates_adaptative"`). Skipped if it already exists.
Output : 
For each cluster:
1. Reads `sequences.a3m` and normalizes sequence names
2. Copies templates from `templates_complete/` to `<cluster_name>/` with cleaned filenames,
   and updates the corresponding MSA sequence names to match.
3. Overwrites `sequences.a3m` with the cleaned names.
4. Writes a `corresponding_name_template.csv` mapping original to cleaned names.
"""
function clean_msa_template_names(clusters,out_folder;cluster_name::String="templates_adaptative")

    for clust in 1:clusters
        corresponding_name=DataFrames.DataFrame(Ref_name=String[],New_name_msa=String[],New_name_template=String[]) # Create an empty DF

        # MSA

        cluster_folder = joinpath(out_folder, "cluster_$(clust)")
        msa_path = joinpath(cluster_folder, "sequences.a3m")

        #1- read the msa
        ids, sequences = read_a3m(msa_path)
        query_id = ids[1]
        query_seq = sequences[1]
        
        #2 - Clean name
        clean_ids=String[]
        for name in ids
            if name == query_id
                push!(clean_ids, query_id)
            elseif startswith(basename(name), "AF")
                push!(clean_ids, uppercase(name))
            else
                pdb_id=first(splitext(basename(name)))
                chain_id = split(basename(name), "_")
                if length(chain_id) > 1
                    chain_id = chain_id[end]
                    
                else
                    chain_id = "a"
                end
                push!(clean_ids, lowercase("$(pdb_id)_$(chain_id)"))
            end
        end
        
        cluster_template_folder = joinpath(cluster_folder, "templates_complete")
        if isdir(joinpath(cluster_folder, cluster_name))
            println("Folder already clean")
            continue
        end
        cluster_template_clean_folder = mkdir(joinpath(cluster_folder, cluster_name))
        
        templates = glob("*",cluster_template_folder)
        isempty(templates) && error("No templates found in $cluster_template_folder")
        i=0
        for template in templates 
            i+=1
            #4- Change template name 
            if startswith(basename(template), "AF")
                clean_template_name = "t00$(i)_a"
                file_name="t00$(i).pdb"
                push!(corresponding_name,(basename(template),clean_template_name,file_name))
                name=first(splitext(basename(template)))
                for i in 1:length(clean_ids)
                    if clean_ids[i]== uppercase(name)
                        clean_ids[i] = clean_template_name
                    end
                end
            else
                clean_template_name = basename(template)
                checks = split(clean_template_name, "_")
                name=split(clean_template_name, "_")[1]
                if length(checks) < 2
                    clean_template_name = name*"_a"
                end
                base_name = replace(name, ".cif" => "")
                file_name=lowercase(base_name)*".cif"
            end
            cp(template, joinpath(cluster_template_clean_folder,file_name); force=true)
        end    

        
        #5- Rewrite MSA with clean names
        open(msa_path, "w") do io
            for (id, seq) in zip(clean_ids, sequences)
                println(io, ">$id")
                println(io, seq)
            end
        end
            
        CSV.write(cluster_folder*"/corresponding_name_template.csv", corresponding_name) 
    end
end