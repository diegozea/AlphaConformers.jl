using AlphaConformers
using Documenter

DocMeta.setdocmeta!(AlphaConformers, :DocTestSetup, :(using AlphaConformers); recursive=true)

makedocs(;
    modules=[AlphaConformers],
    authors="Diego Javier Zea <diegozea@gmail.com> and contributors",
    repo="https://github.com/diegozea/AlphaConformers.jl/blob/{commit}{path}#{line}",
    sitename="AlphaConformers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://diegozea.github.io/AlphaConformers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/diegozea/AlphaConformers.jl",
    devbranch="main",
)
