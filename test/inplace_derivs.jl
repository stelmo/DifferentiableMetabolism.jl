@testset "Inplace differentiable functions using GECKO" begin
    optimizer = Ipopt.Optimizer
    modifications = []

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

    #: Differentiate an optimal solution
    diffmodel = with_parameters(gm, rid_enzyme)
    update_Q!(diffmodel, _ -> spdiagm(fill(0.1, length(diffmodel.var_ids))))

    x_ref, dx_ref = differentiate(diffmodel, Ipopt.Optimizer; scale_output=false)

    rf = ReferenceFunctions(diffmodel)

    # rescale_factors = scaling_factor(rf.E(diffmodel.θ), rf.d(diffmodel.θ));
    rescale_factors = fill(1.0, 10)

    var_derivs, par_derivs = make_inplace_derivatives_with_equality_scaling(rf);

    diffmodel.analytic_var_derivs = (x, ν, λ, θ) -> var_derivs(x, ν, λ, θ, rescale_factors)
    diffmodel.analytic_par_derivs = (x, ν, λ, θ) -> par_derivs(x, ν, λ, θ, rescale_factors)
    
    x = zeros(rf.nx)
    ν = zeros(rf.neq)
    λ = zeros(rf.nineq)
    A = spzeros(rf.nx + rf.neq + rf.nineq, rf.nx + rf.neq + rf.nineq)
    B = spzeros(rf.nx + rf.neq + rf.nineq, rf.nθ)
    dx = zeros(rf.nx, rf.nθ)

    differentiate!(x, ν, λ, A, B, dx, diffmodel, Ipopt.Optimizer)

    # test if reproduceable solutions
    @test all([
        isapprox(x_ref[i], x[i]; atol = TEST_TOLERANCE) for i in eachindex(x)
    ])

    @test all([
        isapprox(dx_ref[i], dx[i]; atol = TEST_TOLERANCE) for i in eachindex(dx)
    ])
end
