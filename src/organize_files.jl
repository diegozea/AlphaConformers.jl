
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

function get_all_predictions(output_dir::String)
    clusters=glob("cluster*",output_dir)
    @show length(clusters)
    dic_pred_struct=Dict()
    for clu in clusters
        predictions=glob("*.pdb",joinpath(clu,"af","predictions","sequences","models"))
        @show length(predictions)
        for pred in predictions
            name=basename(clu)*"_"*basename(pred)
            structure=MIToS.PDB.read_file(pred, MIToS.PDB.PDBFile, group="ATOM")
            dic_pred_struct[name]=structure
        end
    end
    @show length(dic_pred_struct)
    return dic_pred_struct
end  

function compare_struct(dic_pred_struct::Dict, query_struct::String,cutoff_min::Float64,cutoff_max::Float64)
    cluster_close_query=Dict()
    @show "Read query structure"
    if endswith(query_struct, ".cif")
        query_structure=MIToS.PDB.read_file(query_struct, MIToS.PDB.MMCIFFile, group="ATOM")
    else
        query_structure=MIToS.PDB.read_file(query_struct, MIToS.PDB.PDBFile, group="ATOM")
    end
    
    for (name, structure) in dic_pred_struct
        aligned_a, aligned_b, matches, rmsd, coverage, identity=structural_alignment(query_structure, structure)
        @show rmsd
        if rmsd < cutoff_max && rmsd > cutoff_min
            @show name, rmsd
            cluster_close_query[name]=(rmsd)
        end
    end
    @show "Number of predictions close to query: ", length(cluster_close_query)
    @show first(cluster_close_query,5)
    return cluster_close_query
end

function found_uniprot_structure(
    search_results::DataFrames.DataFrame,
    sifts_uniprot_mapping::DataFrames.DataFrame,
    query_pdb_code::String,
    query_chain_code::String)
    @show query_pdb_code, query_chain_code

    query_uniprot = only(get_uniprot_acc(
        sifts_uniprot_mapping,
        query_pdb_code,
        query_chain_code
    ))

    @show query_uniprot

    query_structures = get_pdb_codes(
        sifts_uniprot_mapping,
        String(query_uniprot)
    )

    @show query_structures

    # 👉 colonnes directement
    pdbs   = String.(query_structures.PDB)
    chains = String.(query_structures.CHAIN)

    # 👉 construction vectorisée
    pdb_names = Set(uppercase.(pdbs) .* ".cif_" .* uppercase.(chains))
    pdb_prefixes_upper = Set(uppercase.(pdbs))
    pdb_prefixes_lower = Set(lowercase.(pdbs))

    af_pattern = "AF-$(query_uniprot)-"

    uniprot_result=filter(row -> begin
        t = row.target

        !occursin(af_pattern, t) &&
        !(t in pdb_names) &&
        !any(startswith(t, p) for p in pdb_prefixes_upper) &&
        !any(startswith(t, p) for p in pdb_prefixes_lower)
        end, search_results)

    return uniprot_result
end

function compare_alternative_structures(search_results)
    range_rmsd=[]
    file_analysed=Set{String}()
    for result in eachrow(search_results)
        file_name=result.target
        @show file_name
        if file_name in file_analysed
            continue
        end
        push!(file_analysed, file_name)
        query_pdb_code=String(split(file_name,"_")[1])
        query_chain_code=String(split(file_name,"_")[end])
        center_uniprot=found_uniprot_structure(search_results, sifts_uniprot_mapping, query_pdb_code, query_chain_code)
        @show "Number of alternative structures found for $file_name: ", length(center_uniprot)
        if isempty(center_uniprot)
            continue
        end
        max_rmsd=0.0
        for i in 1:nrow(center_uniprot)
            @show "Comparing $file_name to ", center_uniprot[i, :file]
            push!(file_analysed, center_uniprot[i, :file])
            struct1= MIToS.PDB.read_file(center_uniprot[i, :file], MIToS.PDB.PDBFile, group="ATOM")
            for j in i+1:nrow(center_uniprot)
                struct2=MIToS.PDB.read_file(center_uniprot[j, :file], MIToS.PDB.PDBFile, group="ATOM")

                aligned_a, aligned_b, matches, rmsd, coverage, identity=structural_alignment(struct1, struct2)
                if rmsd > max_rmsd
                    @show "RMSD between ", center_uniprot[i, :file], " and ", center_uniprot[j, :file], " : ", rmsd
                    max_rmsd = rmsd
                end
            end
        end
        push!(range_rmsd, max_rmsd)
    end
    @show "RMSD range between alternative structures: ", range_rmsd
    return minimum(range_rmsd), maximum(range_rmsd)
end

function found_best_prediction(output_dir::String,query_struct::String,sifts_uniprot_mapping)
    dic_pred_struct=get_all_predictions(output_dir)
    cluster_close_query=compare_struct(dic_pred_struct, query_struct,0.0,4.0)

    if isdir(joinpath(output_dir,"target_db_results"))
        @show "Comparing to known structures in target_db_results folder"
        m8file=glob("*.m8",joinpath(output_dir,"target_db_results")) 
        if !isempty(m8file) && isfile(m8file[1])
            @show "Reading foldseek search results from ", m8file[1]
            search_results=read_foldseek_search_results(m8file[1])
            @show first(search_results,5)
            alternative_files=known_uniprot_structures(sifts_uniprot_mapping,search_results)
            @show alternative_files
            full_m8file=glob("*.m8",joinpath(output_dir,"fullpdb_mmcif_files_results")) 
            @show "Reading foldseek all_search results PDB from ", full_m8file[1]
            all_search_results=read_foldseek_search_results(full_m8file[1])
            @show first(all_search_results,5)
            remain_results = similar(search_results, 0)
            for fname in alternative_files
                @show fname
                if startswith(uppercase(fname), uppercase(basename(query_struct)))
                    @show "Skipping query structure ", fname
                    continue
                end
                prot_name=String(split(fname,"_")[1])
                chain_code=String(split(fname,"_")[2])
                filtered_results=filter(r -> startswith(r.target, prot_name) && endswith(r.target, chain_code), all_search_results)
                @show filtered_results
                append!(remain_results, filtered_results)
            end
            if nrow(remain_results) < 2
                @warn "No alternative structures found in foldseek results. Skipping foldseek comparison."
                return cluster_close_query, Dict()
            end
            @show "Comparing alternative structures to query structure"
            @show "Number of alternative structures found: ", nrow(remain_results)
            min_rmsd, max_rmsd = compare_alternative_structures(remain_results)
            cluster_close_objectif=compare_struct(dic_pred_struct, query_struct,min_rmsd,max_rmsd)
            return cluster_close_query, cluster_close_objectif
        else
            @warn "No .m8 file found in target_db_results. Skipping foldseek results."
        end
        
    end
    return cluster_close_query, Dict()
end