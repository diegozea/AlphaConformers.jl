using AlphaConformers
using Test
using TestItemRunner

@testset "AlphaConformers.jl" begin
    include("aqua.jl")
end

@run_package_tests verbose=true
