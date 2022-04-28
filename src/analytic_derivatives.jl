"""
$(TYPEDSIGNATURES)

Creates analytic derivative functions of the `diffmodel` internals, using
symbolic variables. Note, this function can take some time to construct the
derivatives, but substantially speeds up repeated calls to
[`differentiate`](@ref).
"""
function make_derivatives_unscaled(diffmodel::DifferentiableModel)
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

    sparse_F(z) = [
        dm.Q(z[θidxs]) * z[xidxs] + dm.c(z[θidxs]) - dm.E(z[θidxs])' * z[νidxs] -
        dm.M(z[θidxs])' * z[λidxs]
        dm.E(z[θidxs]) * z[xidxs] - dm.d(z[θidxs])
        z[λidxs] .* (dm.M(z[θidxs]) * z[xidxs] - dm.h(z[θidxs]))
    ]
    sz = [sx; sν; sλ; sθ]

    #TODO only get jacobian of theta
    sj = sparse(Symbolics.jacobian(sparse_F(sz), sz)[:, θidxs])
    
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


# struct ReferenceFunctions
#     Q::Function
#     c::Function
#     E::Function
#     d::Function
#     M::Function
#     h::Function
#     nx::Int
#     neq::Int
#     nineq::Int
#     nθ::Int    
# end

# function ReferenceFunctions(dm::DifferentiableModel)
#     ReferenceFunctions(
#         dm.Q,
#         dm.c,
#         dm.E,
#         dm.d,
#         dm.M,
#         dm.h,
#         length(dm.c(dm.θ)),
#         length(dm.d(dm.θ)),
#         length(dm.h(dm.θ)),
#         length(dm.θ),
#     )
# end

# function make_derivatives_scaled(rf::ReferenceFunctions)

#     xidxs = 1:rf.nx
#     νidxs = last(xidxs) .+ (1:rf.neq)
#     λidxs = last(νidxs) .+ (1:rf.nineq)
#     θidxs = last(λidxs) .+ (1:rf.nθ)
#     rfidxs = last(θidxs) .+ (1:rf.neq)

#     Symbolics.@variables begin
#         sx[1:length(xidxs)]
#         sν[1:length(νidxs)]
#         sλ[1:length(λidxs)]
#         sθ[1:length(θidxs)]
#         srowfacts[1:length(νidxs)] #  row rescaler
#     end

#     sparse_F(z) = [
#         rf.Q(z[θidxs]) * z[xidxs] + rf.c(z[θidxs]) - (z[rfidxs] .* rf.E(z[θidxs]))' * z[νidxs] - rf.M(z[θidxs])' * z[λidxs]
#         (z[rfidxs] .* rf.E(z[θidxs])) * z[xidxs] - rf.d(z[θidxs])
#         z[λidxs] .* (rf.M(z[θidxs]) * z[xidxs] - rf.h(z[θidxs]))
#     ]
#     sz = [sx; sν; sλ; sθ; srowfacts]

#     sj = sparse(Symbolics.jacobian(sparse_F(sz), sz)[:, θidxs])
#     (nr, nc) = size(sj)

#     f_expr = build_function(sj, [sx; sν; sλ; sθ; srowfacts])
#     myf = eval(last(f_expr))
#     param_derivs = (out, x, ν, λ, θ, rfs) -> reshape(myf([x; ν; λ; θ; rfs]), nr, nc)

#     var_derivs = (x, ν, λ, θ, rfs) -> [
#         rf.Q(θ) -(rfs .* rf.E(θ))' -rf.M(θ)'
#         (rfs .* rf.E(θ)) spzeros(size(rf.E(θ), 1), length(ν)) spzeros(size(rf.E(θ), 1), length(λ))
#         spdiagm(λ)*rf.M(θ) spzeros(size(rf.M(θ), 1), length(ν)) spdiagm(rf.M(θ) * x - rf.h(θ))
#     ]
    
#     return var_derivs, param_derivs
# end