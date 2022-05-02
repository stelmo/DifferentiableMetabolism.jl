@testset "Differentiable Thermodynamic GECKO" begin
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

    rid_dg0 = Dict("r4" => -1.0, "r5" => -10.0, "r6" => 10.0)

    mid_concentration = Dict(
        "m6" => 0.001,
        "m3" => 0.001,
        "m2" => 1.0,
        "m4" => 0.01,
        "m5" => 0.05,
        "m1" => 1.0,
    )

    #: Solve normal gecko model without thermodynamic effects 
    opt_model = flux_balance_analysis(
        gm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )
    gecko_fluxes = flux_dict(gm, opt_model)
    gecko_gps = gene_product_dict(gm, opt_model)

    #: Differentiate model normal gecko
    diffmodel = with_parameters(
        gm,
        rid_enzyme,
        rid_dg0,
        mid_concentration;
        ignore_reaction_ids = ["r6"],
    )

    x, dx = differentiate(
        diffmodel,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )
    sol = Dict(diffmodel.var_ids .=> x)

    #=
    sanity check
    reaction 4 should require more enzyme due to it being closer to equilibrium
    grr reaction 4: (1 * g2 && 3 * g3)
    =#
    @test sol["g2"] > gecko_gps["g2"]
    @test sol["g3"] > gecko_gps["g3"]

    # check if x, dx can be reproduced
    x_ref = [
        1.2380726767460195
        2.476145353495478
        1.2380726767411319
        1.2380726767489294
        1.2380726767312122
        2.476145353495478
        1.2380726767460195
        0.12380726766577349
        0.061967128585559375
        0.2035881413551585
        0.03537351102455509
    ]
    @test all([isapprox(x_ref[i], x[i]; atol = TEST_TOLERANCE) for i in eachindex(x)])

    dx_ref = [
        0.194554 0.123807 0.681638 -0.0 3.44365e-8 0.341868 -0.341868 -0.341868 -3.44365e-8
        0.194554 0.123807 0.681638 -0.0 3.44365e-8 0.341868 -0.341868 -0.341868 -3.44365e-8
        0.194554 0.123807 0.681638 -0.0 3.44365e-8 0.341868 -0.341868 -0.341868 -3.44365e-8
        0.194554 0.123807 0.681638 -0.0 3.44365e-8 0.341868 -0.341868 -0.341868 -3.44365e-8
        0.194554 0.123807 0.681638 -0.0 3.44365e-8 0.341868 -0.341868 -0.341868 -3.44365e-8
        0.194554 0.123807 0.681638 -0.0 3.44365e-8 0.341868 -0.341868 -0.341868 -3.44365e-8
        0.194554 0.123807 0.681638 0.0 3.44365e-8 0.341868 -0.341868 -0.341868 -3.44365e-8
        0.194554 -0.876193 0.681638 -0.0 3.44365e-8 0.341868 -0.341868 -0.341868 -3.44365e-8
        0.194554 0.123807 -0.318362 -0.0 3.44365e-8 -0.159671 0.159671 0.159671 -3.44365e-8
        0.107679 0.123807 -0.231486 -0.0 1.90594e-8 -0.116099 0.116099 0.116099 -1.90594e-8
        -0.805446 0.123807 0.681638 -0.0 -1.42566e-7 0.341868 -0.341868 -0.341868 1.42566e-7
    ]
    @test all([
        isapprox(dx_ref[i], dx[1:11, :][i]; atol = TEST_TOLERANCE) for
        i in eachindex(dx_ref)
    ])

    # test if automatic and symbolic derivatives are the same
    make_derivatives(diffmodel)
    _, dx_sym = differentiate(diffmodel, Tulip.Optimizer; use_analytic_nonmutating = true)
    @test all([
        isapprox(dx_sym[i], dx[i]; atol = TEST_TOLERANCE) for i in eachindex(dx_sym)
    ])
end
