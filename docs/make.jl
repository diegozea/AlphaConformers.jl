using Pkg

# The docs manifest is ignored, so make the documented package available when
# the docs command is run from a fresh checkout.
cd(@__DIR__) do
    Pkg.develop(PackageSpec(path = ".."))
    Pkg.instantiate()
end

using AlphaConformers
using Documenter

DocMeta.setdocmeta!(
    AlphaConformers,
    :DocTestSetup,
    :(using AlphaConformers);
    recursive = true,
)

makedocs(;
    modules = [AlphaConformers],
    authors = "Diego Javier Zea <diegozea@gmail.com> and contributors",
    repo = Documenter.Remotes.GitHub("diegozea", "AlphaConformers.jl"),
    sitename = "AlphaConformers.jl",
    workdir = joinpath(@__DIR__, ".."),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://diegozea.github.io/AlphaConformers.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Database Setup" => "databases.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/diegozea/AlphaConformers.jl.git", devbranch = "main")
