@testset "Differentiable SMOMENT" begin
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

    #: Get SMoment solution for comparison
    fluxes = flux_balance_analysis_dict(
        smm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    #: Diffeentiate model 
    diffmodel = with_parameters(smm, rid_enzyme)

    x, dx = differentiate(
        diffmodel,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    # test if COBREXA Smoment is the same as the differentiable one 
    sol = Dict(reactions(smm) .=> x)
    @test isapprox(sol["r1"], fluxes["r1"]; atol = TEST_TOLERANCE)
    @test isapprox(sol["r6"], fluxes["r6"]; atol = TEST_TOLERANCE)

    # test if automatic and symbolic derivatives are the same
    make_derivatives(diffmodel)
    _, dx_sym = differentiate(diffmodel, Tulip.Optimizer; use_analytic = true)
    @test all([
        isapprox(dx_sym[i], dx[i]; atol = TEST_TOLERANCE) for i in eachindex(dx_sym)
    ])

    # test if reference solution is attained
    dx_ref = [
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
        0.251908 0.160305 0.587786
    ]
    @test all([
        isapprox(dx_ref[i], dx[1:7, :][i]; atol = TEST_TOLERANCE) for i in eachindex(dx_ref)
    ])
end
