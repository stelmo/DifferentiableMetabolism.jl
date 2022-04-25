@testset "Differentiable GECKO" begin
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

    #: Get classic GECKO solution
    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )
    gecko_fluxes = flux_dict(gm, opt_model)
    gecko_gps = gene_product_dict(gm, opt_model)
    gene_product_mass_group_dict(gm, opt_model)

    #: Differentiate an optimal solution
    diffmodel = with_parameters(gm, rid_enzyme; scale_equality = false)

    x_noscaling, dx_noscaling = differentiate(
        diffmodel,
        Tulip.Optimizer;
        use_analytic = false,
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    (slb, sub) = check_scaling(diffmodel)

    @test isapprox(slb, 1.845098040014257; atol = TEST_TOLERANCE)
    @test isapprox(sub, 2.146128035678238; atol = TEST_TOLERANCE)

    diffmodel = with_parameters(gm, rid_enzyme; scale_equality = true)

    x_scaling, dx_scaling = differentiate(
        diffmodel,
        Tulip.Optimizer;
        use_analytic = false,
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test all([
        isapprox(dx_noscaling[i], dx_scaling[i]; atol = TEST_TOLERANCE) for
        i in eachindex(dx_scaling)
    ])

    @test all([
        isapprox(x_noscaling[i], x_scaling[i]; atol = TEST_TOLERANCE) for
        i in eachindex(x_scaling)
    ])
end