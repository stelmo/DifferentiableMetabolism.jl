"""
$(TYPEDSIGNATURES)

Solve and differentiate an optimization problem using the optimality conditions.
The output can be scaled relative to the parameters and the solved variables
with `scale_output`. *Optimizer* modifications (from COBREXA.jl) can be supplied
through `modifications`. Analytic derivatives of the optimality conditions can
be used by setting `use_analytic` to true. 

Internally calls [`differentiate!`](@ref), the in-place variant of this function, 
constructs the arguments internally.
"""
function differentiate(
    diffmodel::DifferentiableModel,
    optimizer;
    use_analytic = false,
    scale_output = true,
    modifications = [],
)

    x = zeros(length(diffmodel.var_ids))
    ν = zeros(length(diffmodel.d(diffmodel.θ)))
    λ = zeros(length(diffmodel.h(diffmodel.θ)))
    nA = length(x) + length(ν) + length(λ)
    A = spzeros(nA, nA)
    B = zeros(nA, length(diffmodel.param_ids))
    dx = zeros(size(A, 1), size(B, 2))

    differentiate!(
        x,
        ν,
        λ,
        A,
        B,
        dx,
        diffmodel,
        optimizer;
        modifications,
        use_analytic,
        scale_output,
    )

    return x, dx
end

"""
$(TYPEDSIGNATURES)

Differentiates `diffmodel` using in-place functions as much as possible. The meaning of 
the variables are:
- `x`: primal variables
- `ν`: equality dual variables
- `λ`: inequality dual variables
- `A`: variable derivatives
- `B`: parameter derivitives
- `dx`: the sensitivities of the variables
- `diffmodel`: the model to differantiate
- `optimizer`: the solver used in the forward pass
- `modifications`: COBREXA functions forwarded to the *solver*
- `use_analytic`: use the analytic derivatives (need to be supplied)
- `scale_output`: flag if the primal variable sensitivities should be scaled by `dx/dy => dlog(x)/dlog(y)`
"""
function differentiate!(
    x,
    ν,
    λ,
    A,
    B,
    dx,
    diffmodel::DifferentiableModel,
    optimizer;
    modifications = [],
    use_analytic = false,
    scale_output = true,
)
    DifferentiableMetabolism._solve_model!(x, ν, λ, diffmodel, optimizer; modifications)

    DifferentiableMetabolism._differentiate_kkt!(
        x,
        ν,
        λ,
        A,
        B,
        dx,
        diffmodel;
        use_analytic,
    )

    #: Scale dx/dy => dlog(x)/dlog(y)
    scale_output && _scale_derivatives(x, dx, diffmodel)

    nothing
end

"""
$(TYPEDSIGNATURES)

Solve the optimization problem, aka the forward pass.

This function is called by [`differentiate!`](@ref). 
"""
function _solve_model!(
    x,
    ν,
    λ,
    diffmodel::DifferentiableModel,
    optimizer;
    modifications = [],
)
    #: forward pass, solve the optimization problem
    opt_model = JuMP.Model(optimizer)
    set_silent(opt_model)
    @variable(opt_model, z[1:length(diffmodel.var_ids)])

    if all(diffmodel.Q(diffmodel.θ) .== 0) # is LP
        @objective(opt_model, Min, diffmodel.c(diffmodel.θ)' * z)
    else
        @objective(
            opt_model,
            Min,
            0.5 * z' * diffmodel.Q(diffmodel.θ) * z + diffmodel.c(diffmodel.θ)' * z
        )
    end

    @constraint(opt_model, eq, diffmodel.E(diffmodel.θ) * z .== diffmodel.d(diffmodel.θ))
    @constraint(opt_model, ineq, diffmodel.M(diffmodel.θ) * z .<= diffmodel.h(diffmodel.θ))

    # apply the modifications
    for mod in modifications
        mod(nothing, opt_model)
    end

    optimize!(opt_model)

    if termination_status(opt_model) ∉ [JuMP.OPTIMAL, JuMP.LOCALLY_SOLVED]
        throw(DomainError(termination_status(opt_model), " model not solved optimally!"))
    end

    #: differentiate the optimal solution
    x .= value.(opt_model[:z])
    ν .= dual.(opt_model[:eq])
    λ .= dual.(opt_model[:ineq])

    nothing
end

"""
$(TYPEDSIGNATURES)

Implicitly differentiate a convex quadratic or linear program using the KKT
conditions. Suppose `F(z(θ), θ) = 0` represents the optimality (KKT) conditions
of some optimization problem specified by `diffmodel`, where `z` is a vector of
the primal, `x`, and dual, `ν` and `λ`, solutions. Then, mutate `A = ∂₁F(z(θ),
θ)`, `B = ∂₂F(z(θ), θ)`, and the optimal solution.

This function is called by [`differentiate!`](@ref). 
"""
function _differentiate_kkt!(
    x,
    ν,
    λ,
    A,
    B,
    dx,
    diffmodel::DifferentiableModel;
    use_analytic = false,
)
    #=
    The analytic approaches are much faster than using automatic
    differentiation, but more labor intensive because the derivative of the
    parameters with respect to the KKT function needs to be manually supplied.
    Ideally, use symbolic code to generate the derivatives for you.

    However, the autodiff works for any model and the user does not have to
    supply any analytic derivatives. However, any sparse arrays encoded
    in the model structure must be cast to dense arrays for ForwardDiff to
    work. The densification restriction can be dropped once #481 here 
    https://github.com/JuliaDiff/ForwardDiff.jl/pull/481 gets merged into 
    ForwardDiff.
    =#
    if use_analytic
        A .= diffmodel.analytic_var_derivs(x, ν, λ, diffmodel.θ)
        B .= diffmodel.analytic_par_derivs(x, ν, λ, diffmodel.θ)
    else
        xidxs = 1:length(x)
        νidxs = last(xidxs) .+ (1:length(ν))
        λidxs = last(νidxs) .+ (1:length(λ))
        θidxs = last(λidxs) .+ (1:length(diffmodel.θ))

        kkt(z) = [
            Array(diffmodel.Q(z[θidxs])) * z[xidxs] + Array(diffmodel.c(z[θidxs])) -
            Array(diffmodel.E(z[θidxs]))' * z[νidxs] -
            Array(diffmodel.M(z[θidxs]))' * z[λidxs]
            Array(diffmodel.E(z[θidxs])) * z[xidxs] - Array(diffmodel.d(z[θidxs]))
            diagm(z[λidxs]) *
            (Array(diffmodel.M(z[θidxs])) * z[xidxs] - Array(diffmodel.h(z[θidxs])))
        ]

        J = ForwardDiff.jacobian(kkt, [x; ν; λ; diffmodel.θ])
        A .= J[:, 1:last(λidxs)] #TODO this can be analytic by default
        B .= J[:, (1+last(λidxs)):end]
    end

    #= 
    Solve the linear system of equations.
    This is the most time intensive operation, accounting for 75% 
    of the run time on iML1515. Iterative solvers, and other indirect 
    approaches did not speed up the solving process. Separating out 
    the factorization also does not help that much, because the pattern 
    changes, requiring frequent refactorizations. Sadly, Pardiso is no 
    longer free for academics. Hence, the standard solve is used. 
    =# 
    _linear_solver(A, B, dx)

    nothing
end

"""
$(TYPEDSIGNATURES)

Separate linear solve to have the ability to swap it our for something faster.
"""
function _linear_solver(A, B, dx)
    ldiv!(dx, lu(-A), B)
    nothing
end

"""
$(TYPEDSIGNATURES)

Scale the derivatives: `dx/dy => dlog(x)/dlog(y)`. Note, only scales 
the primal variables, *not* the duals.
"""
function _scale_derivatives(x, dx, diffmodel::DifferentiableModel)
    for i = 1:length(x)
        for j = 1:size(dx, 2)
            dx[i, j] *= diffmodel.θ[j] / x[i]
        end
    end
    nothing
end
