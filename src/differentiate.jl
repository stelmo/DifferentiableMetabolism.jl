"""
    differentiate(
        diffmodel::DifferentiableModel,
        optimizer;
        use_analytic = false,
        scale_equality = false,
        scale_output = true,
        modifications = [],
        regularizer = 0.0,
    )

Solve and differentiate an optimization problem using the optimality conditions.
Optionally, scale the equality constraints with `scale_equality`. The output can
be scaled relative to the parameters and the solved variables with
`scale_output`. *Optimizer* modifications (from COBREXA.jl) can be supplied
through `modifications`. Analytic derivatives of the optimality conditions can
be used by setting `use_analytic` to true.

Internally calls [`_differentiate_kkt`](@ref).
"""
function differentiate(
    diffmodel::DifferentiableModel,
    optimizer;
    use_analytic = false,
    scale_equality = false,
    scale_output = true,
    modifications = [],
    regularizer = 0.0,
)

    A, B, x = _differentiate_kkt(
        diffmodel,
        optimizer;
        modifications,
        scale_equality,
        use_analytic,
        regularizer,
    )

    dx = -sparse(A) \ Array(B) # no method for sparse \ sparse
    dx = dx[1:length(diffmodel.var_ids), :] # only return derivatives of variables, not the duals

    #: Scale dx/dy => dlog(x)/dlog(y)
    if scale_output
        scaled_dx = similar(dx)
        for i = 1:size(dx, 1)
            for j = 1:size(dx, 2)
                scaled_dx[i, j] = diffmodel.θ[j] / x[i] * dx[i, j]
            end
        end
        return x, scaled_dx
    else
        return x, dx
    end
end

"""
    _differentiate_kkt(
        diffmodel::DifferentiableModel,
        optimizer;
        modifications = [],
        use_analytic = false,
        scale_equality = false,
        regularizer = 0.0,
    )

Implicitly differentiate a convex quadratic or linear program using the KKT
conditions. Suppose `F(z(θ), θ) = 0` represents the optimality (KKT) conditions
of some optimization problem specified by `diffmodel`, where `z` is a vector of
the primal and dual solutions. Then, return `A = ∂₁F(z(θ), θ)`, `B = ∂₂F(z(θ),
θ)`, and the optimal solution.

This function is called by [`differentiate`](@ref). 
"""
function _differentiate_kkt(
    diffmodel::DifferentiableModel,
    optimizer;
    modifications = [],
    use_analytic = false,
    scale_equality = false,
    regularizer = 0.0,
)
    #: forward pass, solve the optimization problem
    E = diffmodel.E
    h = diffmodel.h
    M = diffmodel.M
    d = diffmodel.d
    c = diffmodel.c
    Q = diffmodel.Q
    θ = diffmodel.θ
    regQ = spdiagm(fill(regularizer, length(diffmodel.var_ids)))

    if scale_equality
        row_factors = scaling_factor(E(θ), d(θ))
    else
        row_factors = fill(1.0, size(E(θ), 1))
    end

    opt_model = Model(optimizer)
    set_silent(opt_model)
    @variable(opt_model, x[1:size(E(θ), 2)])

    if all(Q(θ) + regQ .== 0)
        @objective(opt_model, Min, c(θ)' * x)
    else
        @objective(opt_model, Min, 0.5 * x' * (Q(θ) + regQ) * x + c(θ)' * x)
    end

    @constraint(opt_model, eq, row_factors .* E(θ) * x .== row_factors .* d(θ))
    @constraint(opt_model, ineq, M(θ) * x .<= h(θ))

    # apply the modifications
    for mod in modifications
        mod(nothing, opt_model)
    end

    optimize!(opt_model)

    @assert(termination_status(opt_model) in [JuMP.OPTIMAL, JuMP.LOCALLY_SOLVED])

    #: differentiate the optimal solution
    x = value.(opt_model[:x])
    ν = dual.(opt_model[:eq])
    λ = dual.(opt_model[:ineq])

    if use_analytic
        #=
        This is much faster than using automatic differentiation, but more labor 
        intensive because the derivative of the parameters with respect to the KKT function 
        needs to be manually supplied.
        =#
        A = [
            Q(θ)+regQ -(row_factors .* E(θ))' -M(θ)'
            row_factors.*E(θ) zeros(size(E(θ), 1), length(ν)) zeros(size(E(θ), 1), length(λ))
            diagm(λ)*M(θ) zeros(size(M(θ), 1), length(ν)) diagm(M(θ) * x - h(θ))
        ]
        B = diffmodel.param_derivs(x, ν, λ, θ, regularizer)
    else
        #=
        This is slower but works for any model and the user does not have to supply any analytic derivatives.
        However, the any sparse arrays encoded in the model structure must be cast to dense arrays for ForwardDiff to work. 
        =#
        F(x, ν, λ, θ) = [
            (Array(Q(θ)) + Array(regQ)) * x + Array(c(θ)) -
            (row_factors .* Array(E(θ)))' * ν - Array(M(θ))' * λ
            (row_factors .* Array(E(θ))) * x - (row_factors .* Array(d(θ)))
            diagm(λ) * (Array(M(θ)) * x - Array(h(θ)))
        ]

        A = [ForwardDiff.jacobian(x -> F(x, ν, λ, θ), x) ForwardDiff.jacobian(
            ν -> F(x, ν, λ, θ),
            ν,
        ) ForwardDiff.jacobian(λ -> F(x, ν, λ, θ), λ)]
        B = ForwardDiff.jacobian(θ -> F(x, ν, λ, θ), θ)
    end

    return A, B, x
end
