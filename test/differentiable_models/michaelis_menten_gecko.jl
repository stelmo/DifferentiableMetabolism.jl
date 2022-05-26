@testset "Differentiable Michaelis Menten GECKO" begin
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

    rid_km = Dict(
        "r4" => Dict("m3" => 0.007, "m4" => 0.02, "m5" => 0.03),
        "r5" => Dict("m2" => 0.4, "m4" => 0.02, "m6" => 0.01),
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
        rid_enzyme;
        rid_dg0,
        rid_km,
        mid_concentration,
        ignore_reaction_ids = ["r6"],
    )

    x, dx = differentiate(
        diffmodel,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )
    sol = Dict(diffmodel.var_ids .=> x)

    # check if x, dx can be reproduced
    x_ref = [
        0.09530050777800092
        0.19060101555594944
        0.09530050777775771
        0.09530050777819793
        0.09530050777821945
        0.19060101555594944
        0.09530050777800092
        0.009530050776567447
        0.050561076614846864
        0.1911648757697115
        0.07896329185025613
    ]
    @test all([isapprox(x_ref[i], x[i]; atol = TEST_TOLERANCE) for i in eachindex(x)])

    dx_ref = [
        0.434298 0.00953005 0.556172 -0.0 0.119806 0.719682 0.644083 0.224761 0.419322
        0.434298 0.00953005 0.556172 -0.0 0.119806 0.719682 0.644083 0.224761 0.419322
        0.434298 0.00953005 0.556172 -0.0 0.119806 0.719682 0.644083 0.224761 0.419322
        0.434298 0.00953005 0.556172 -0.0 0.119806 0.719682 0.644083 0.224761 0.419322
        0.434298 0.00953005 0.556172 -0.0 0.119806 0.719682 0.644083 0.224761 0.419322
        0.434298 0.00953005 0.556172 -0.0 0.119806 0.719682 0.644083 0.224761 0.419322
        0.434298 0.00953005 0.556172 0.0 0.119806 0.719682 0.644083 0.224761 0.419322
        0.434298 -0.99047 0.556172 -0.0 0.119806 0.719682 0.644083 0.224761 0.419322
        0.434298 0.00953005 -0.443828 -0.0 0.119806 -0.57431 0.239962 -0.179361 0.419322
        0.227766 0.00953005 -0.237296 -0.0 0.0628321 -0.307059 0.124016 -0.0958966 0.219912
        -0.565702 0.00953005 0.556172 -0.0 -0.156056 0.719682 -0.321434 0.224761 -0.546195
    ]
    @test all([
        isapprox(dx_ref[i], dx[1:11, :][i]; atol = TEST_TOLERANCE) for
        i in eachindex(dx_ref)
    ])

    # test if automatic and symbolic derivatives are the same
    make_derivatives(diffmodel)
    _, dx_sym = differentiate(
        diffmodel,
        Tulip.Optimizer;
        use_analytic = true,
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test all([
        isapprox(dx_sym[i], dx[i]; atol = TEST_TOLERANCE) for i in eachindex(dx_ref)
    ])
end
