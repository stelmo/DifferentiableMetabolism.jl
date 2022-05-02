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
    diffmodel = with_parameters(gm, rid_enzyme;)

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

    # test if automatic and symbolic derivatives are the same
    make_derivatives(diffmodel)
    x_anal, dx_anal = differentiate(diffmodel, Tulip.Optimizer; use_analytic = true)
    @test all([
        isapprox(dx_auto[i], dx_anal[i]; atol = TEST_TOLERANCE) for i in eachindex(dx_anal)
    ])

    #: Add a regularizer and test QP 
    update_Q!(diffmodel, _ -> spdiagm(fill(0.1, length(diffmodel.var_ids))))
    x_qp, dx_qp = differentiate(diffmodel, Ipopt.Optimizer; use_analytic = false)

    # test if reproduceable solutions
    x_qp_ref = [
        0.7677550153758821
        1.5355100307517642
        0.7677550153758821
        0.7677550153758821
        0.7677550153758821
        1.5355100307517642
        0.7677550153758821
        0.07677550153758822
        0.025591833845862735
        0.08774343032867224
        0.02193585758216806
    ]
    @test all([
        isapprox(x_qp_ref[i], x_qp[i]; atol = TEST_TOLERANCE) for i in eachindex(x_qp)
    ])

    dx_qp_ref = [
        0.000376045 0.00153551 0.00192548
        0.000376045 0.00153551 0.00192548
        0.000376045 0.00153551 0.00192548
        0.000376045 0.00153551 0.00192548
        0.000376045 0.00153551 0.00192548
        0.000376045 0.00153551 0.00192548
        0.000376045 0.00153551 0.00192548
        0.000376045 -0.998464 0.00192548
        0.000376045 0.00153551 -0.998075
        -0.124624 0.00153551 -0.873075
        -0.999624 0.00153551 0.00192548
    ]
    @test all([
        isapprox(dx_qp_ref[i], dx_qp[1:11, :][i]; atol = TEST_TOLERANCE) for i in eachindex(dx_qp_ref)
    ])
end
