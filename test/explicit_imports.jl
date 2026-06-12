@testitem "explicit imports lint" begin
    using AlphaConformers
    import ExplicitImports

    ExplicitImports.test_explicit_imports(
        AlphaConformers;
        all_qualified_accesses_are_public=false,
    )
end
