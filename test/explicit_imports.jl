@testitem "explicit imports lint" begin
    using AlphaConformers
    import ExplicitImports

    # `all_qualified_accesses_*` are relaxed because the package intentionally reaches re-exported
    # symbols through their umbrella package (e.g. `Plots.mm`, whose owner is `Measures`).
    ExplicitImports.test_explicit_imports(
        AlphaConformers;
        all_qualified_accesses_are_public = false,
        all_qualified_accesses_via_owners = false,
    )
end
