"""
    differentiate(
        diffmodel::DifferentiableModel,
        optimizer;
        use_analytic = false,
        scale_output = true,
        modifications = [],
    )

Solve and differentiate an optimization problem using the optimality conditions.
The output can be scaled relative to the parameters and the solved variables
with `scale_output`. *Optimizer* modifications (from COBREXA.jl) can be supplied
through `modifications`. Analytic derivatives of the optimality conditions can
be used by setting `use_analytic` to true.

Internally calls [`_differentiate_kkt`](@ref).
"""
function differentiate(
    diffmodel::DifferentiableModel,
    optimizer;
    use_analytic = false,
    scale_output = true,
    modifications = [],
)

    A, B, x = _differentiate_kkt(diffmodel, optimizer; modifications, use_analytic)

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
)
    #: forward pass, solve the optimization problem
    Q = diffmodel.Q
    c = diffmodel.c
    E = diffmodel.E
    d = diffmodel.d
    M = diffmodel.M
    h = diffmodel.h
    θ = diffmodel.θ

    opt_model = Model(optimizer)
    set_silent(opt_model)
    @variable(opt_model, x[1:size(E(θ), 2)])

    if all(Q(θ) .== 0)
        @objective(opt_model, Min, c(θ)' * x)
    else
        @objective(opt_model, Min, 0.5 * x' * Q(θ) * x + c(θ)' * x)
    end

    @constraint(opt_model, eq, E(θ) * x .== d(θ))
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
        intensive because the derivative of the parameters with respect to the
        KKT function needs to be manually supplied.
        =#
        A = [
            Q(θ) -E(θ)' -M(θ)'
            E(θ) zeros(size(E(θ), 1), length(ν)) zeros(size(E(θ), 1), length(λ))
            diagm(λ)*M(θ) zeros(size(M(θ), 1), length(ν)) diagm(M(θ) * x - h(θ))
        ]
        B = diffmodel.param_derivs(x, ν, λ, θ)
    else
        #=
        Note, if the internal structures of diffmodel are changes, then
        auto_derivs MUST be updated. 
        =#
        A, B = diffmodel.auto_derivs(x, ν, λ, θ)
    end

    return A, B, x
end
