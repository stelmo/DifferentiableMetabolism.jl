@testset "Gecko" begin
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
    diffmodel = with_parameters(
        gm,
        rid_enzyme;
        analytic_parameter_derivatives = derivative_of_enzyme_equality(gm, rid_enzyme),
    )

    x_auto, dx_auto = differentiate(
        diffmodel,
        Tulip.Optimizer;
        use_analytic = false,
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    # test if solution is the same between gecko and differentiable gecko
    sol = Dict(diffmodel.var_ids .=> x_auto)
    @test isapprox(sol["g1"], gecko_gps["g1"]; atol = TEST_TOLERANCE)
    @test isapprox(sol["r3#forward#1"], gecko_fluxes["r3"]; atol = TEST_TOLERANCE)

    # test if automatic and manual derivatives are the same
    x_anal, dx_anal = differentiate(diffmodel, optimizer; use_analytic = true)
    @test all([
        isapprox(dx_auto[i], dx_anal[i]; atol = TEST_TOLERANCE) for i in eachindex(dx_anal)
    ])

    #: Add a regularizer and test QP
    x_qp, dx_qp =
        differentiate(diffmodel, Ipopt.Optimizer; use_analytic = false, regularizer = 0.1)

    # test if reproduceable solutions
    x_qp_ref = [
        1.6030534486107917
        3.2061068972215834
        1.6030534486107917
        1.6030534486107917
        1.6030534486107917
        3.2061068972215834
        1.6030534486107917
        0.16030534486107917
        0.05343511495369306
        0.1832061084126619
        0.04580152710316548
    ]
    @test all([
        isapprox(x_qp_ref[i], x_qp[i]; atol = TEST_TOLERANCE) for i in eachindex(x_qp)
    ])

    dx_qp_ref = [
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 -0.839695 0.587786
        0.251908 0.160305 -0.412214
        0.126908 0.160305 -0.287214
        -0.748092 0.160305 0.587786
    ]
    @test all([
        isapprox(dx_qp_ref[i], dx_qp[i]; atol = TEST_TOLERANCE) for i in eachindex(dx_qp)
    ])
end
