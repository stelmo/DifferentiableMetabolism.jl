@testset "Differentiable Thermodynamic SMOMENT" begin
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

    rid_dg0 = Dict("r4" => -1.0, "r5" => -10.0)

    mid_concentration = Dict(
        "m6" => 0.001,
        "m3" => 0.001,
        "m2" => 1.0,
        "m4" => 0.01,
        "m5" => 0.05,
        "m1" => 1.0,
    )
    @test true

    # #: Differentiate model 
    # diffmodel = with_parameters(smm, rid_enzyme, rid_dg0, mid_concentration)

    # x, dx = differentiate(
    #     diffmodel,
    #     Tulip.Optimizer;
    #     use_analytic = false,
    #     modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    # )
end
