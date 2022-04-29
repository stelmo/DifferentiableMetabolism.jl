"""
$(TYPEDSIGNATURES)

Differentiates `diffmodel` using in-place functions as much as possible.
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
)
    _differentiate_kkt!(x, ν, λ, A, B, diffmodel, optimizer; modifications)

    # lua = lu(-A)
    # ldiv!(dx, lua, B)
    # n = length(diffmodel.var_ids)
    # for (i, b) in enumerate(eachcol(B))   
    #     dx[:, i] .= (lua\b)[1:n]
    # end
    dx .= (-sparse(A) \ Array(B))[1:length(diffmodel.var_ids), :]

    nothing
end

"""
$(TYPEDSIGNATURES)

Implicit differentiation, in-place.
"""
function _differentiate_kkt!(
    x,
    ν,
    λ,
    A,
    B,
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
        @objective(opt_model, Min, 0.5 * z' * diffmodel.Q(diffmodel.θ) * z + diffmodel.c(diffmodel.θ)' * z)
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

    A .= diffmodel.analytic_var_derivs(x, ν, λ, diffmodel.θ)
    B .= diffmodel.analytic_par_derivs(x, ν, λ, diffmodel.θ)

    nothing
end