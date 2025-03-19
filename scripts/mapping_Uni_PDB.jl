#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l ncpus=4
#PBS -l mem=16g
#PBS -l host=node48
#PBS -l walltime=300:00:00
#PBS -j oe
=#

import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")

using AlphaConformers
using BioStructures
using DataFrames
using MIToS.PDB
using MIToS.SIFTS
using MIToS.MSA
using JSON3
using Statistics
import CSV
import MIToS
using OrderedCollections
using DataStructures
using Plots
import BioAlignments
using Clustering 
using Statistics
using LinearAlgebra
using PairwiseListMatrices
using StatsBase

#Lie les 3 csv ensembles
"""
Create one file with all of the information form PFAM CATH and UNIPROT from SIFT databases 

Take in input the three file il a DataFrame format and output one DataFrame
"""
function join_information(sift_pfam_mapping::DataFrame,sift_cath_mapping::DataFrame,sift_uniprot_mapping::DataFrame)

    # 1 Extraire RES_BEG et RES_END depuis sift_uniprot_mapping
    uniprot_res = select(sift_uniprot_mapping, [:PDB, :CHAIN, :RES_BEG, :RES_END, :SP_PRIMARY])
    # 2 Compter le nombre de PFAM_ID par (PDB, CHAIN, SP_PRIMARY)
    pfam_counts = combine(groupby(sift_pfam_mapping, [:PDB, :CHAIN, :SP_PRIMARY])) do subdf
        (PFAM_NB = size(subdf, 1),)
    end
    # 3 Compter le nombre de CATH_ID par (PDB, CHAIN, SP_PRIMARY)
    cath_counts = combine(groupby(sift_cath_mapping, [:PDB, :CHAIN, :SP_PRIMARY])) do subdf
        (CATH_NB = size(subdf, 1),)
    end
    # 4 Joindre les deux tables pour créer le DataFrame final
    sift_join_file = outerjoin(uniprot_res, pfam_counts, on=[:PDB, :CHAIN, :SP_PRIMARY])
    sift_join_file = outerjoin(sift_join_file, cath_counts, on=[:PDB, :CHAIN, :SP_PRIMARY],)
    # 5 Remplacer les valeurs `missing` par 0 (au cas où un PDB n'a pas de PFAM ou CATH)
    sift_join_file.PFAM_NB .= coalesce.(sift_join_file.PFAM_NB, 0)
    sift_join_file.CATH_NB .= coalesce.(sift_join_file.CATH_NB, 0)
    # 6 Trier par SP_PRIMARY
    sift_join_file = sort(sift_join_file, :SP_PRIMARY)

    return sift_join_file 
end 

#Recupere information du header utile
"""
Use the function getpdbdescription to have the details about the pdb file 

Take in input a PDB code in lowercase and output the information in a vecteur 
We will have the code pdb, the release date, the method of extraction and the resolution if we have the information 
(can be missing if the method is NMR)
"""
function get_header_pdb(code_pdb::String)
    header_info=getpdbdescription(code_pdb) # Retourne un DICT
    #Extraire les informations
    ## Récupérer la date de publication
    date_info_p = get(header_info, "rcsb_accession_info", Dict())
    date_info_d = get(date_info_p, "initial_release_date", nothing)
    if date_info_d !== nothing
        date_info_t = split(date_info_d, "T")
        date_info = date_info_t[1]
    else
        return nothing
    end
    ## Récupérer la méthode
    exp_info = get(header_info, "rcsb_entry_info", Dict())
    method_info = get(exp_info, "experimental_method", nothing)
    ## Récupérer la résolution
    res_info_p = get(exp_info, "resolution_combined", nothing)
    #Prepare the output
    if res_info_p !== nothing && !isempty(res_info_p)
        res_info = res_info_p[1]
        return (code_pdb, date_info, method_info, res_info)
    else 
        return (code_pdb, date_info, method_info,missing)
    end
end

#Enregistre information dans un DF
"""
From the Dataframe with all the information of Pfam CATH and Uniprot we had the information form the header 
We use the function get_header_pdb

Take in input the Dataframe from join_information and output a Dataframe with the header information added 
"""
function get_pdb_information(sift_join_file::DataFrame)
    #Creer DF vide
    df_pdb_reso_resi = DataFrame(PDB = String[], CHAIN = String[], RES_BEG=Union{Int64,Missing}[], RES_END=Union{Int64,Missing}[],PFAM_NB=Int64[],CATH_NB=Int64[], RESOLUTION = Union{Float64, Missing}[], METHOD = String[],DATE_RELEASE=SubString{String}[],UNIPROT=String[]) # Prendre en compte si resolution manquante
    i=0   
    #Pour chaque fichier pdb      
    for row in eachrow(sift_join_file)  
        # Recupere code uniprot associé 
        uni_acc=row.SP_PRIMARY
        #Recuperer information du header 
        code_pdb, date_info, method_info, res_info = get_header_pdb(String(row.PDB))
        #Remplis le DF
        push!(df_pdb_reso_resi,(code_pdb,String(row.CHAIN),row.RES_BEG,row.RES_END,row.PFAM_NB,row.CATH_NB,res_info,method_info,date_info,String(uni_acc))) 
    end 
    return df_pdb_reso_res
end

#Recuperer les ligands 
"""
Get the type of ligand that can bound with each pdb, if no ligand we put missing

Take in input the DataFrame from BioLip and the DataFrame created with get_pdb_information and add to it in output 
the ligand (We are looking the type not the quantities) 
"""
function get_ligand_information(df_biolip_5_first_columns::DataFrame,df_pdb_reso_resi::DataFrame)
    # Recupere une ligne pour chaque pdb --> group les ligands 
    df_ligands = combine(groupby(df_biolip_5_first_columns, [:PDB, :CHAIN]), 
                     :LIGAND => (x -> join(unique(x), ", ")) => :LIGANDS)
    # Fusionner avec le DF principal
    df_final = leftjoin(df_pdb_reso_resi, df_ligands, on=[:PDB, :CHAIN])
    # Remplacer les valeurs manquantes par "missing"
    df_final.LIGANDS .= coalesce.(df_final.LIGANDS, missing)
    return df_final
end 

#Recuperer le lien entre les indices Uniprot et PDB
"""
Do the correspondance between Uniprot and PDB index. Use SIFTS mapping so we need to download the .xml file.
We handle the error if the file couldn't be install 

Take in input a row for the DataFrame get_ligand_information and the path to a folder to save temporary file
output a dictionnaire with the mapping Uniprot => PDB
"""
function get_sift_mapping(row::DataFrameRow,folder_temporary_path::String)
    pdb=row.PDB
    #Chemin où est enregistrer le fichier xml 
    siftsfile=joinpath(folder_temporary_path,row.PDB*".xml.gz")
    if !isfile(siftsfile) #Verifier que fichier n'est pas déja télécharger 
        println("The file doesn't exist. Downloading the file XML...")
        siftsfile = downloadsifts(uppercase(pdb),filename=joinpath(folder_temporary_path,row.PDB*".xml.gz" ))
        #Si le fichier n'a pas pu etre télécharger 
        pdbcode=row.PDB
        if siftsfile == nothing 
            try 
                sleep(60)
                filename=joinpath(folder_temporary_path,lowercase(pdbcode)*".xml.gz" )
                siftsfile=download_file(string("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/", 
                lowercase(pdbcode), ".xml.gz"), filename)
            catch e
                println("❌ Error when downloading $pdbcode with 2nd link: ", e)
            end
            if siftsfile == nothing
                return nothing
            end
            println("✅ Found $pdbcode with the 2nd link")
        else 
            println("✅ Found $pdbcode ")  
        end
    end
    #Do the mapping 
    siftsmap = siftsmapping(  # Retourne un Dictionnaire avec coordonnée Uniprot => cordonnée PDB 
        siftsfile,
        dbUniProt,
        String(row.UNIPROT),
        dbPDB,
        String(pdb), # SIFTS stores PDB identifiers in lowercase
        chain = String(row.CHAIN), # In this example we are only using the chain A of the PDB
        missings = false,
    ) # Residues without coordinates aren't used in the mapping
    return siftsmap
end

#Faire le lien entre les indices Uniprots et les résidues
"""
Make the link between the Uniprot index and the residues. Use the fonction get_sift_mapping and MIToS.PDB.downloadpdb

Take in input one row form the DataFrame get_ligand_information and the path to a folder to save temporary file
output a dictionnaire with the mapping Uniprot => Residues
"""
function get_uniprot_mapping_residues(row::DataFrameRow,folder_temporary_path::String)
    # Recupere le mapping entre Uniprot => PDB
    mapping=get_sift_mapping(row,folder_temporary_path) #OrderedDict Uniprot => PDB
    if mapping == nothing
        return nothing # Si le fichier n'a pas pu etre télécharger 
    end
    # Recupere les residues PDB 
    pdbfile=joinpath(folder_temporary_path,uppercase(row.PDB)*".pdb.gz" )
    if !isfile(pdbfile) #Verifier que fichier n'est pas déja télécharger 
        println("The file doesn't exist. Downloading the PDB ...")
        pdb=row.PDB
        try 
            pdbfile = MIToS.PDB.downloadpdb(uppercase(row.PDB), format = PDBFile,filename=joinpath(folder_temporary_path,uppercase(row.PDB)*".pdb.gz" )) # Sinon ne reconnais pas le package --> weird   
        catch e 
            println("❌ Error when downloading the PDB file $pdb",e)
            return nothing
        end   
    end
    residues_1ivo = read(pdbfile, PDBFile)
    chain_residues_pdb = MIToS.PDB.residuesdict(residues_1ivo, "1", String(row.CHAIN), "ATOM") #Retourne dictionnaire avec position => residues details 
    #Pour tous les residues trouver l'index Uniprot
    residues_uniprot = OrderedDict()
    for id_uni in keys(mapping)
        try 
            id_pdb=mapping[string(id_uni)] #Doit etre un string pour les comparer 
            residue=chain_residues_pdb[id_pdb]
            residues_uniprot[string(id_uni)]=residue # Output Ordred Dict Uniprot => Residue 
        catch e 
            pdb=row.PDB
            println("❌ Key $id_uni not found for $pdb")
        end 
    end
    return residues_uniprot
end

"""
Compare the residues form PDB and Uniprot to identify the number of mutation and the missing residues

Take in input the DataFrame from get_ligand_information and the path to a folder to save temporary file
Output the Dataframe with two new colonnes MUTATION and MISSING_RESIDUES
"""
function  count_mutation_pdb_uni(df_uniprot::DataFrame,folder_temporary_path::String)
    #Ajoute colonne avec mutation 
    df_uniprot.MUTATION .= 0 
    #Ajoute colonne avec residue manquant
    df_uniprot.MISSING_RESIDUES .= 0
    # Pour chaque fichier pdb 
    for row in eachrow(df_uniprot)
        mut=0
        pdb=row.PDB
        # Recupere le fichier SIFT
        siftsfile = downloadsifts(uppercase(pdb),filename=joinpath(folder_temporary_path,row.PDB*".xml.gz" ))
        pdbcode=row.PDB
        if siftsfile == nothing
            try 
                sleep(60) #Sleep to make sure we don't try to many time 
                filename=joinpath(folder_temporary_path,lowercase(pdbcode)*".xml.gz" )
                siftsfile=download_file(string("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/", 
                lowercase(pdbcode), ".xml.gz"), filename)
            catch e
                println("❌ Error when downloading $pdbcode with the 2nd link: ", e)
            end
            if siftsfile == nothing
                return df_uniprot #Pour les fichiers qui n'ont pas réussi à etre installé 
            end
            println("✅ Found $pdbcode with 2nd link")
        else 
            println("✅ Found $pdbcode ")  
        end

        residue_data = read(siftsfile, SIFTSXML)
        #Boucle pour chaque residues
        for i in 1:length(residue_data)
            #Recupere les informations
            residues=read(siftsfile, SIFTSXML)[i]
            res_uni=get(residues, dbUniProt, :name, "missing")
            id_uni=get(residues, dbUniProt, :id, "missing")
            res_pdb=get(residues, dbPDB, :name, "")
            chain_pdb=get(residues, dbPDB, :chain, "")
            #Verifie que l'on est sur la bonne chain et bon Uniprot
            #Si pas dans uniprot on ignore
            if row.CHAIN==chain_pdb && row.UNIPROT==id_uni
                res_pdb_n=three2residue(res_pdb) # Retourne un MIToS.MSA.Residue
                if res_uni != string(res_pdb_n) # Si pas même residue alors mutation
                    mut=mut+1
                end
            end  
        end
        #Ajoute valeur dans df 
        row.MUTATION = mut
        #Compte le nombre de missing residues 
        sifts_1ivo = read(siftsfile, SIFTSXML, chain = String(row.CHAIN))
        mis=count([res.missing for res in sifts_1ivo]) # Retourne un vecteur booleen donc compte le nombre de 1
        row.MISSING_RESIDUES=mis
    end
    return df_uniprot
end

"""
Compare for each uniprot accession number the pdb file associate with structure alignment 

Take in input the DataFrame from count_mutation_pdb_uni and the path to a folder to save temporary file
Output a Vector with the RMSD between each PDB --> its the half of the Correlation matrix
"""
function calculate_RMSD_uniprot(df_uniprot::DataFrame,folder_temporary_path::String)
    files = df_uniprot.PDB .* "_" .* df_uniprot.CHAIN  # Liste des fichiers
    n = length(files)  # Nombre de fichiers
    rmsd_list = []
    coverage_list[]
    for i in 1:n
        for j in (i+1):n  # Évite les doublons 
            #file1,file2=files[i],files[j]
            dict1=get_uniprot_mapping_residues(df_uniprot[i,:],folder_temporary_path)
            dict2 = get_uniprot_mapping_residues(df_uniprot[j,:],folder_temporary_path)
            if dict1 == nothing  || dict2 == nothing
                rmsd = missing
            else
                common_keys = Set(keys(dict1)) ∩ Set(keys(dict2))
                push!(coverage_list,common_keys)
                if length(common_keys)==0
                    rmsd=missing
                else 

                    positions = SortedSet(parse.(Int,common_keys))  # Récupère les positions communes (même clé dans tous les fichiers)
                    # Extraction des résidus pour chaque fichier sous forme de liste ordonnée
                    residues1 = [dict1[string(pos)] for pos in positions] 
                    residues2 = [dict2[string(pos)] for pos in positions]
                    try  
                        # Superposition et calcul du RMSD
                        superimposed_1, superimposed_2,rmsd = superimpose(residues1, residues2) # Le 2e élément est la RMSD
                    catch e 
                        println(common_keys)
                        println("❌ Erreur lors de la superposition des structures ", e)
                        rmsd = missing
                    end 
                end
            end
            #p=plot(superimposed_1, label="1", alpha=0.5)
            #plot!(superimposed_2, label="2", alpha=0.5)
            #savefig(p,file1*file2*"2.png")
            push!(rmsd_list, rmsd)  # Stocke la valeur
        end
    end
    plm=PairwiseListMatrix(coverage_list,false,100.0)
    nplm=setlabels(plm,files)
    println(nplm)
    return rmsd_list,files
end

"""
Create the cluster for different CutOff : 1.5 2 and 4

Take in input the Vector with the RMSD from calculate_RMSD_uniprot
Output a Dataframe with the cluster number for each CutOff
"""
function cah(rmsd_list,files)
    #Creer la matrice symétrique
    plm=PairwiseListMatrix(rmsd_list,false,0.0)
    nplm=setlabels(plm,files)
    # Supprimer les lignes et colonnes qui contiennent encore des `missing`
    rows_with_missing = findall(row -> any(ismissing, row), eachrow(nplm))
    df_rmsd_2 = nplm[setdiff(1:end, rows_with_missing), setdiff(1:end, rows_with_missing)]

    mat = Matrix{Float64}(df_rmsd_2) # Conversion en Matrice
    clustering_result = hclust(mat;linkage=:average)#Faire le clustering

    #Découper les cluster en fonction de different cutoff
    cutoffs = [1.5, 2.0, 4.0] # Liste des différents cutoffs à tester

    #Creation du DF output
    pdb_values = [split(file, "_")[1] for file in files]
    chain_values = [split(file, "_")[2] for file in files]
    final_df = DataFrame("PDB" => pdb_values,"CHAIN"=>chain_values)

    for cutoff in cutoffs  # Découper les clusters à différents niveaux de cutoff
        cluster_assignments = cutree(clustering_result, h=cutoff)
        #Enregistre les résultats
        df_clusters = DataFrame()
        file_i = files[setdiff(1:end, rows_with_missing)]
        df_clusters.PDB = [split(f, "_")[1] for f in file_i]  
        df_clusters.CHAIN = [split(f, "_")[2] for f in file_i] 
        df_clusters[!, Symbol("Cluster_$cutoff")] = cluster_assignments  # Ajout des clusters

        # Faire un `left join` pour inclure les PDB qui n'ont pas pu etre analyser 
        final_df = leftjoin(final_df, df_clusters, on=on=[:PDB,:CHAIN])
    end
    return final_df
end

"""
For each uniprot accession we are adding information about the number of mutation and the missing residues, we calculate the RMSD 
between every pdb form the pdb accesion to compare there stucture and than cluster them in cluster base on there similarity

Take in input the DataFrame from get_ligand_information
Output the DataFrame with new columns MUTATION MISSING_RESIDUES Cluster_1.5 Cluster_2.0 Cluster_4.0
"""
function clustering_each_uniprot_acc(df_final::DataFrame)
    df_completed=DataFrame(UNIPROT=String[],MUTATION=Int64[]) #Pour pouvoir checker lors de la premiere itération
    i=0
    #Regarde pour chaque uniprot acc 
    for row in eachrow(df_final)
        uniprot=row.UNIPROT
        println(uniprot)
        is_present = true
        if i!=0 # Regarde pas pour la premiere itération 
            is_present= uniprot in df_completed.UNIPROT # Regarde si Uniprot est déja dans le DF pour ne pas faire plusieurs fois même calcul 
        end
        #if !is_present  || i==0 # Si uniprot pas present dans le DF où si c'est la premiere itération 
        if uniprot =="A0A024B7W1"    
            #Creation d'un dossier vide
            df_completed=mktempdir() do temp_dir
                df_uniprot = filter(row -> row.UNIPROT == uniprot, df_final)# Pour ce centrer par uniprot 
                #println(df_uniprot)

                # Verifier les differences entre Uniprot et PDB
                df_uniprot=count_mutation_pdb_uni(df_uniprot,temp_dir)
                #println(first(df_uniprot,10))

                #Calcul RMSD entre chaque fichier pdb 
                rmsd_list,files=calculate_RMSD_uniprot(df_uniprot,temp_dir)
                #println(df_rmsd)

                #Creer les clusters
                df_cluster=cah(rmsd_list,files)
                #println(df_cluster)

                #Regroupe les resultats dans un DF
                select!(df_uniprot, Not(["RES_END", "RES_BEG"]))
                df_result = leftjoin(df_uniprot, df_cluster, on=[:PDB,:CHAIN])
                println(df_result)
                if i ==0
                    df_completed=df_result
                else 
                    df_completed=vcat(df_completed, df_result)
                end
                #println(first(df_completed,10))
                return df_completed
            end
        end
        i=i+1
    end
    return df_completed
end

"""
Function to analyse the different cutoff from th cluster and chose the more releatable one 

Take in input the DataFrame from clustering_each_uniprot_acc
Output a DataFrame to know for each Uniprot accession if we have the apo and holo form and for with cluster there where different
"""
function check_apo_holo_cluster(df_completed::DataFrame)
    check_cluster=DataFrame(UNIPROT=String[],HOLO_APO=Bool[],BEST_CUTOFF=Union{Float64,Missing}[])
    for row in eachrow(df_completed)
        uniprot=row.UNIPROT
        #println(uniprot)
        #Verifie qu'on a pas déja fait 
        is_present= uniprot in check_cluster.UNIPROT
        if !is_present
            df_uniprot = select(filter(row -> row.UNIPROT == uniprot, df_completed), ["PDB", "LIGANDS", "Cluster_1.5", "Cluster_2.0","Cluster_4.0"])# Pour ce centrer par uniprot et garder que information qui nous interesse 
            #println(df_uniprot)
            is_apo = any(ismissing, df_uniprot.LIGANDS)
            if is_apo 
                # Filtrer les fichiers avec et sans ligands
                with_ligands = filter(row -> !ismissing(row.LIGANDS), df_uniprot)
                without_ligands = filter(row -> ismissing(row.LIGANDS), df_uniprot)
                result = 9
                # Comparer les clusters pour chaque cutoff
                for cutoff in [1.5, 2.0, 4.0]
                    # Récupérer les clusters pour chaque cutoff
                    cluster_col = Symbol("Cluster_$(cutoff)")
                    #println(cluster_col)
                    # Comparer les fichiers sans ligands et avec ligands
                    for file_with_ligands in eachrow(with_ligands)
                        for file_without_ligands in eachrow(without_ligands)
                            if !ismissing(file_with_ligands[cluster_col]) && 
                                !ismissing(file_without_ligands[cluster_col]) && 
                                file_with_ligands[cluster_col] != file_without_ligands[cluster_col]
                                # Si les fichiers sont dans des clusters différents, prendre la valeur la plus petite
                                if  cutoff<result
                                    result = cutoff
                                end
                            end
                        end
                    end
                end
                push!(check_cluster,(uniprot,true,result))
                
            else 
                push!(check_cluster,(uniprot,false,missing))
            end 
        end
    end
    return check_cluster
end
### MAIN ####

println("START ! ")

## Get Uniprot 
sift_uniprot_mapping=get_uniprot_mapping()
## Get pfam
sifts_file_path_pfam="pdb_chain_pfam.csv"
sift_pfam_mapping=DataFrames.DataFrame(CSV.File(sifts_file_path_pfam,
        comment="#", missingstring=["", "None"])) #Output DF with PDB CHAIN SP_PRIMARY PFAM_ID COVERAGE

##Get CATH 
sift_file_path_cath="pdb_chain_cath_uniprot.csv"
sift_cath_mapping=DataFrames.DataFrame(CSV.File(sift_file_path_cath,
        comment="#", missingstring=["", "None"])) #Output DF with PDB CHAIN SP_PRIMARY CATH_ID

"""
## Get biolip
file_path_biolip="BioLiP.csv"
df_biolip=CSV.File(file_path_biolip, delim='\t') |> DataFrame

# Garder uniquement les 5 premières colonnes
df_biolip_5_first_columns = select(df_biolip, 1:5)
rename!(df_biolip_5_first_columns, :Column2=> "CHAIN")
rename!(df_biolip_5_first_columns, :Column3=> "RESOLUTION")
rename!(df_biolip_5_first_columns, :Column4=> "SITE")
rename!(df_biolip_5_first_columns, :Column5=> "LIGAND")

# Sauvegarder en CSV
CSV.write("Biolip_5_first_columns.csv", df_biolip_5_first_columns)
println(first(df_biolip_5_first_columns,20))
println(size(df_biolip_5_first_columns))

## Get 5 first columns of biolip
file_path_biolip="Biolip_5_first_columns.csv"
df_biolip_5_first_columns=DataFrames.DataFrame(CSV.File(file_path_biolip,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(df_biolip_5_first_columns,20))
println(size(df_biolip_5_first_columns))
rename!(df_biolip_5_first_columns, Symbol.(strip.(string.(names(df_biolip_5_first_columns))))) #Enlever les espaces dans le nom des colonnes 


## Join PFAM and CATH
sift_join_file=join_information(sift_pfam_mapping,sift_cath_mapping,sift_uniprot_mapping) #Output DF with PDB CHAIN RES_BEG RES_END SP_PRIMARY PFAM_NB CATH_NB
println(first(sift_join_file,20))
println(size(sift_join_file))

## Recuperer information pdb 
df_pdb_reso_resi=get_pdb_information(sift_join_file) # Output DF with PDB CHAIN RES_BEG RES_END SP_PRIMARY PFAM_NB CATH_NB
println(first(df_pdb_reso_resi,20))
println(size(df_pdb_reso_resi))


##Recuperer information sur le ligand 
df_final=get_ligand_information(df_biolip_5_first_columns,df_pdb_reso_resi)
println(first(df_final,50))
println(size(df_final))

CSV.write("pdb_information_details.csv", df_final)
println("END !") 

"""
# get back the final df 
file_path_df_final="pdb_information_details_1000.csv"
df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND
println(first(df_final,20))
println(size(df_final))


#Clustering entre même proteine 
df_completed=clustering_each_uniprot_acc(df_final)
println(first(df_completed,20))
println(size(df_completed))

CSV.write("pdb_information_details_final_1000.csv", df_completed)
#Verifie que holo et apo form sont dans des clusters differents
check_cluster=check_apo_holo_cluster(df_completed)
println(first(check_cluster,20))
value_counts = countmap(check_cluster.BEST_CUTOFF)  # Compte les occurrences
most_frequent_value = argmax(value_counts) 
println("Le cutoff qui separe le mieux la forme apo et holo est : ", most_frequent_value)
frequency = count(x -> !ismissing(x) && x == most_frequent_value, check_cluster.BEST_CUTOFF)
println("Elle apparaît ", frequency, " fois.")
println(size(check_cluster))
println("Elle separe correctement : ",(frequency/size(check_cluster)[1])*100)

println("END !")