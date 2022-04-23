@testset "Update" begin
    #: Set problem up 
    stdmodel,
    reaction_isozymes,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound = create_test_model()

    gm = make_gecko_model(
        stdmodel;
        reaction_isozymes,
        gene_product_bounds,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )

    rid_enzyme = Dict(
        k => isozyme_to_enzyme(first(v), gene_product_molar_mass; direction = :forward)
        for (k, v) in reaction_isozymes
    )

    #: Differentiate an optimal solution
    diffmodel = with_parameters(gm, rid_enzyme;)

    #: Test if update works
    n_vars = length(diffmodel.c(diffmodel.θ))
    n_eqs = size(diffmodel.E(diffmodel.θ), 1)
    n_ineqs = size(diffmodel.M(diffmodel.θ), 1)
    θ = diffmodel.θ

    _Q = rand(n_vars, n_vars)
    _c = rand(n_vars)
    _E = rand(n_eqs, n_vars)
    _d = rand(n_eqs)
    _M = rand(n_ineqs, n_vars)
    _h = rand(n_ineqs)

    update_Q!(diffmodel, _ -> _Q)
    update_c!(diffmodel, _ -> _c)
    update_E!(diffmodel, _ -> _E)
    update_d!(diffmodel, _ -> _d)
    update_M!(diffmodel, _ -> _M)
    update_h!(diffmodel, _ -> _h)

    @test all(diffmodel.Q(θ) .== _Q)
    @test all(diffmodel.c(θ) .== _c)
    @test all(diffmodel.E(θ) .== _E)
    @test all(diffmodel.d(θ) .== _d)
    @test all(diffmodel.M(θ) .== _M)
    @test all(diffmodel.h(θ) .== _h)
end
