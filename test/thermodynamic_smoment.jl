@testset "Differentiable Thermodynamic SMOMENT" begin
    #: Set problem up 
    stdmodel, reaction_isozymes, _, gene_product_molar_mass, _ = create_test_model()

    smm = make_smoment_model(
        stdmodel;
        reaction_isozyme = Dict(k => first(v) for (k, v) in reaction_isozymes),
        gene_product_molar_mass,
        total_enzyme_capacity = 1.0,
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

    #: Differentiate model 
    diffmodel = with_parameters(
        smm,
        rid_enzyme,
        rid_dg0,
        mid_concentration;
        ignore_reaction_ids = ["r6"],
    )

    x, dx = differentiate(
        diffmodel,
        Tulip.Optimizer;
        use_analytic = false,
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    # test if can reproduce answers
    x_ref = [
        1.2380726841130605
        2.4761453682260397
        1.2380726841129737
        1.2380726841134497
        1.2380726841130556
        2.4761453682260397
        1.2380726841130605
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
    ]
    @test all([isapprox(dx_ref[i], dx[i]; atol = TEST_TOLERANCE) for i in eachindex(dx)])

    # test if automatic and symbolic derivatives are the same
    make_derivatives(diffmodel)
    _, dx_sym = differentiate(diffmodel, Tulip.Optimizer; use_analytic = true)
    @test all([
        isapprox(dx_sym[i], dx[i]; atol = TEST_TOLERANCE) for i in eachindex(dx_sym)
    ])
end
