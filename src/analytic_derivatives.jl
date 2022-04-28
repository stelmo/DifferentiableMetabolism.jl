"""
$(TYPEDSIGNATURES)

Creates analytic derivative functions of the `diffmodel` internals, using
symbolic variables. Note, this function can take some time to construct the
derivatives, but substantially speeds up repeated calls to
[`differentiate`](@ref).
"""
function make_derivatives(diffmodel::DifferentiableModel)
    diffmodel.analytic_par_derivs = make_symbolic_param_derivative(diffmodel)
    diffmodel.analytic_var_derivs = make_analytic_var_derivative(diffmodel)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Creates the analytic derivatives of the parameters using symbolic programming.
"""
function make_symbolic_param_derivative(dm::DifferentiableModel)
    xidxs = 1:length(dm.c(dm.θ))
    νidxs = last(xidxs) .+ (1:length(dm.d(dm.θ)))
    λidxs = last(νidxs) .+ (1:length(dm.h(dm.θ)))
    θidxs = last(λidxs) .+ (1:length(dm.θ))

    Symbolics.@variables begin
        sx[1:length(xidxs)]
        sν[1:length(νidxs)]
        sλ[1:length(λidxs)]
        sθ[1:length(θidxs)]
    end

    sparse_F(z) = Array([
        dm.Q(z[θidxs]) * z[xidxs] + dm.c(z[θidxs]) - dm.E(z[θidxs])' * z[νidxs] -
        dm.M(z[θidxs])' * z[λidxs]
        dm.E(z[θidxs]) * z[xidxs] - dm.d(z[θidxs])
        z[λidxs] .* (dm.M(z[θidxs]) * z[xidxs] - dm.h(z[θidxs]))
    ])
    sz = [sx; sν; sλ; sθ]

    #TODO only get jacobian of theta (slower, see #580 in Symbolics.jl)
    sj = Symbolics.sparsejacobian(sparse_F(sz), sz)[:, θidxs]
    
    f_expr = build_function(sj, [sx; sν; sλ; sθ])
    myf = eval(first(f_expr))

    (x, ν, λ, θ) -> myf([x; ν; λ; θ])
end

"""
$(TYPEDSIGNATURES)

Creates the analytic derivatives of the variables.
"""
make_analytic_var_derivative(dm::DifferentiableModel) =
    (x, ν, λ, θ) -> [
        dm.Q(θ) -dm.E(θ)' -dm.M(θ)'
        dm.E(θ) spzeros(size(dm.E(θ), 1), length(ν)) spzeros(size(dm.E(θ), 1), length(λ))
        spdiagm(λ)*dm.M(θ) spzeros(size(dm.M(θ), 1), length(ν)) spdiagm(dm.M(θ) * x - dm.h(θ))
    ]
