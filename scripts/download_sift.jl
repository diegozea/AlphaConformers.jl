#!/store/EQUIPES/AMIG/MEMBERS/diego.zea/bin/julia19

#=
#PBS -l host=node48
#PBS -l walltime=300:00:00
#PBS -j oe
=#
import Pkg
Pkg.activate("/home/julie.daniel/.julia/environments/v1.11")
cd("/store/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformers.jl/scripts/")

using BioStructures
using MIToS.PDB
using MIToS.SIFTS
using DataFrames
import CSV

output_directory = abspath("/alpha/database/pdb", "sift_files")
isdir(output_directory) || mkdir(output_directory)

output_directory_pdb = abspath("/alpha/database/pdb", "pdb_files")

#get the pdb we want 
file_path_df_final="pdb_information_details.csv"
df_final=DataFrames.DataFrame(CSV.File(file_path_df_final,
comment="#", missingstring=["", "None"])) # Output DF with PDB CHAIN RESOLUTION SITE LIGAND

@info "START"
#For each pdb 
sift_not_download=[]
pdb_not_download=[]
for row in eachrow(df_final)
    pdb=row.PDB
    #Check if xml file exist 
    siftsfile=joinpath(output_directory,pdb*".xml.gz")
    if !isfile(siftsfile) #if it is not already downloaded
        @info "The file doesn't exist. Downloading the file XML..."
        siftsfile = downloadsifts(uppercase(pdb),filename=joinpath(output_directory,pdb*".xml.gz" ))
        #If the file could not be downloaded
        if siftsfile == nothing 
            try #Try with an other link
                sleep(60)
                filename=joinpath(output_directory,lowercase(pdb)*".xml.gz" )
                siftsfile=download_file(string("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/", 
                lowercase(pdb), ".xml.gz"), filename)
            catch e
                @warn "❌ Error when downloading $pdb with 2nd link: "
            end
            if siftsfile == nothing # If xml still couldn't be download rigth it in a txt
                push!(sift_not_download,pdb)
            end
            @info "✅ Found $pdb with the 2nd link"
        else 
            @info "✅ Found $pdb " 
        end
    end

    #Check if the pdb is download
    pdbfile=joinpath(output_directory_pdb,pdb*".pdb")
    if !isfile(pdbfile) #if it is not already downloaded
        push!(pdb_not_download,pdb)
    end
end

#Create a file with the pdb error 
open("Error_downloading_sift_pdb.txt", "w") do io
    println(io, "SIFT file that couldn't be download by download_sift.jl :")
    for item in sift_not_download
        println(io, item)
    end

    println(io, "")  

    println(io, "PDB file that aren't download in pdb_files :")
    for item in pdb_not_download
        println(io, item)
    end
end
@info "END !"