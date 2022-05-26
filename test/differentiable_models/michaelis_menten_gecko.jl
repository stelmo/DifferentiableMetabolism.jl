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
        0.1259560041366511
        0.2519120082733372
        0.12595600413638888
        0.1259560041369218
        0.1259560041365431
        0.2519120082733372
        0.1259560041366511
        0.012595600412036653
        0.08720892830551037
        0.26418189287384747
        0.0051102159284440365
    ]
    @test all([isapprox(x_ref[i], x[i]; atol = TEST_TOLERANCE) for i in eachindex(x)])

    dx_ref = [
        0.0281062   0.0125956   0.959298   -0.0   0.0083131    1.37108    -0.886044   -0.885648   -0.000395867
        0.0281062   0.0125956   0.959298   -0.0   0.0083131    1.37108    -0.886044   -0.885648   -0.000395867
        0.0281062   0.0125956   0.959298   -0.0   0.0083131    1.37108    -0.886044   -0.885648   -0.000395867
        0.0281062   0.0125956   0.959298   -0.0   0.0083131    1.37108    -0.886044   -0.885648   -0.000395867
        0.0281062   0.0125956   0.959298   -0.0   0.0083131    1.37108    -0.886044   -0.885648   -0.000395867
        0.0281062   0.0125956   0.959298   -0.0   0.0083131    1.37108    -0.886044   -0.885648   -0.000395867
        0.0281062   0.0125956   0.959298    0.0   0.0083131    1.37108    -0.886044   -0.885648   -0.000395867
        0.0281062  -0.987404    0.959298   -0.0   0.0083131    1.37108    -0.886044   -0.885648   -0.000395867
        0.0281062   0.0125956  -0.0407018  -0.0   0.0083131   -0.058173    0.0371811   0.0375769  -0.000395867
        0.0184344   0.0125956  -0.03103    -0.0   0.00545244  -0.0443496   0.028388    0.0286477  -0.000259643
       -0.971894    0.0125956   0.959298   -0.0  -0.287462     1.37108    -0.87196    -0.885648    0.0136888
    ]

    @test all([
        isapprox(dx_ref[i], dx[1:11, :][i]; atol = TEST_TOLERANCE_RELAXED) for
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
