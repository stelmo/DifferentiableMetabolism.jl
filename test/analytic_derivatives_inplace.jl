@testset "Inplace differentiable functions using GECKO" begin
    #=
    Note, the non-in-place analytic derivatives are tested 
    in the models. The inplace derivatives get special treatment,
    because they are more of a pain to setup. 
    =#

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

    rf = ReferenceFunctions(diffmodel)

    #: Run once to get the shape of the in place variables A, B
    x = zeros(rf.nx)
    ν = zeros(rf.neq)
    λ = zeros(rf.nineq)
    _A = spzeros(rf.nx + rf.neq + rf.nineq, rf.nx + rf.neq + rf.nineq)
    _B = zeros(rf.nx + rf.neq + rf.nineq, rf.nθ)
    dx = zeros(rf.nx + rf.neq + rf.nineq, rf.nθ)

    make_derivatives(diffmodel)

    differentiate!(
        x,
        ν,
        λ,
        _A,
        nothing,
        _B,
        dx,
        diffmodel, 
        Ipopt.Optimizer; 
        scale_output=false,
        use_analytic_nonmutating = true,
    )

    x_ref = deepcopy(x)
    dx_ref = deepcopy(dx)

    #: Now use inplace variants for all functions
    rescale_factors = scaling_factor(rf.E(diffmodel.θ), rf.d(diffmodel.θ))
    update_E!(diffmodel, θ ->  rescale_factors .* rf.E(θ))
    update_d!(diffmodel, θ ->  rescale_factors .* rf.d(θ))
   
    var_derivs, par_derivs = make_inplace_derivatives_with_equality_scaling(rf);

    diffmodel.analytic_var_derivs = (A, x, ν, λ, θ) -> var_derivs(A, x, ν, λ, θ, rescale_factors)
    diffmodel.analytic_par_derivs = (B, x, ν, λ, θ) -> par_derivs(B, x, ν, λ, θ, rescale_factors)
    
    A = sparse(_A)
    B = sparse(_B)
    fA = lu(_A)

    differentiate!(
        x, 
        ν, 
        λ, 
        A, 
        fA,
        B, 
        dx, 
        diffmodel, 
        Ipopt.Optimizer;
        scale_output=false,
        use_analytic_mutating = true,
    )

    # test if reproduceable solutions
    @test all([
        isapprox(x_ref[i], x[i]; atol = TEST_TOLERANCE) for i in eachindex(x)
    ])

    @test all([
        isapprox(dx_ref[1:11,:][i], dx[1:11,:][i]; atol = TEST_TOLERANCE) for i in eachindex(dx[1:11,:])
    ])

    #: now test if truly in place
    fA_prev_L = deepcopy(fA.L)
    fA_prev_U = deepcopy(fA.U)
    
    diffmodel.θ .*= 1.1
    differentiate!(
        x, 
        ν, 
        λ, 
        A, 
        fA,
        B, 
        dx, 
        diffmodel, 
        Ipopt.Optimizer;
        scale_output=false,
        use_analytic_mutating = true,
    )
    @test !all(fA.L .== fA_prev_L)
    @test !all(fA.U .== fA_prev_U)
end
