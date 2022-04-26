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

    #: Differentiate an optimal solution
    diffmodel = with_parameters(gm, rid_enzyme; scale_equality = false)

    x_noscaling, dx_noscaling = differentiate(
        diffmodel,
        Tulip.Optimizer;
        use_analytic = false,
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    (slb, sub), _, _ = check_scaling(diffmodel)

    @test isapprox(slb, 1.845098040014257; atol = TEST_TOLERANCE)
    @test isapprox(sub, 2.146128035678238; atol = TEST_TOLERANCE)

    #: scale equality
    diffmodel = with_parameters(gm, rid_enzyme; scale_equality = true)

    x_scaling_eq, dx_scaling_eq = differentiate(
        diffmodel,
        Tulip.Optimizer;
        use_analytic = false,
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test all([
        isapprox(dx_noscaling[i], dx_scaling_eq[i]; atol = TEST_TOLERANCE) for
        i in eachindex(dx_scaling)
    ])

    @test all([
        isapprox(x_noscaling[i], x_scaling_eq[i]; atol = TEST_TOLERANCE) for
        i in eachindex(x_scaling)
    ])

    #: scale both Inequality and equality
    diffmodel = with_parameters(gm, rid_enzyme; scale_equality = true, scale_inequality = true)

    x_scaling_eqineq, dx_scaling_eqineq = differentiate(
        diffmodel,
        Tulip.Optimizer;
        use_analytic = false,
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test all([
        isapprox(dx_noscaling[i], dx_scaling_eqineq[i]; atol = TEST_TOLERANCE) for
        i in eachindex(dx_scaling)
    ])

    @test all([
        isapprox(x_noscaling[i], x_scaling_eqineq[i]; atol = TEST_TOLERANCE) for
        i in eachindex(x_scaling)
    ])
end
