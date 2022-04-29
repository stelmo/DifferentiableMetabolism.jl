"""
$(TYPEDEF)

A struct used to store information the original optimization problem.

$(TYPEDFIELDS)
"""
struct ReferenceFunctions
    Q::Function
    c::Function
    E::Function
    d::Function
    M::Function
    h::Function
    nx::Int
    neq::Int
    nineq::Int
    nθ::Int    
end

function Base.show(io::IO, ::MIME"text/plain", m::ReferenceFunctions)
    println(io, "", "A reference function collection.")
end

"""
$(TYPEDSIGNATURES)

A constructor for [`ReferenceFunctions`](@ref).
"""
function ReferenceFunctions(dm::DifferentiableModel)
    ReferenceFunctions(
        dm.Q,
        dm.c,
        dm.E,
        dm.d,
        dm.M,
        dm.h,
        length(dm.c(dm.θ)),
        length(dm.d(dm.θ)),
        length(dm.h(dm.θ)),
        length(dm.θ),
    )
end

"""
$(TYPEDSIGNATURES)

Creates in place derivatives.
"""
function make_inplace_derivatives_with_equality_scaling(rf::ReferenceFunctions)

    param_derivs = make_inplace_param_deriv_with_scaling(rf)
    var_derivs = make_inplace_var_deriv_with_scaling(rf)
    
    return var_derivs, param_derivs
end

function make_inplace_param_deriv_with_scaling(rf::ReferenceFunctions)
    xidxs = 1:rf.nx
    νidxs = last(xidxs) .+ (1:rf.neq)
    λidxs = last(νidxs) .+ (1:rf.nineq)
    θidxs = last(λidxs) .+ (1:rf.nθ)
    rfidxs = last(θidxs) .+ (1:rf.neq)

    Symbolics.@variables begin
        sx[1:length(xidxs)]
        sν[1:length(νidxs)]
        sλ[1:length(λidxs)]
        sθ[1:length(θidxs)]
        srowfacts[1:length(νidxs)] #  row rescaler
    end

    sparse_F(z) = Array([
        rf.Q(z[θidxs]) * z[xidxs] + rf.c(z[θidxs]) - (z[rfidxs] .* rf.E(z[θidxs]))' * z[νidxs] - rf.M(z[θidxs])' * z[λidxs]
        (z[rfidxs] .* rf.E(z[θidxs])) * z[xidxs] - rf.d(z[θidxs])
        z[λidxs] .* (rf.M(z[θidxs]) * z[xidxs] - rf.h(z[θidxs]))
    ])
    sz = [sx; sν; sλ; sθ; srowfacts]
    # sparse_F(sz)

    sj = Symbolics.sparsejacobian(sparse_F(sz), sz)[:, θidxs]

    f_expr = build_function(sj, [sx; sν; sλ; sθ; srowfacts], expression=Val{false})
    myf = first(f_expr)

    (x, ν, λ, θ, rfs) -> reshape(myf([x; ν; λ; θ; rfs]), size(sj)...) #TODO fix this
end

function make_inplace_var_deriv_with_scaling(rf::ReferenceFunctions)
    (x, ν, λ, θ, rfs) -> [
        rf.Q(θ) -(rfs .* rf.E(θ))' -rf.M(θ)'
        (rfs .* rf.E(θ)) spzeros(size(rf.E(θ), 1), length(ν)) spzeros(size(rf.E(θ), 1), length(λ))
        spdiagm(λ)*rf.M(θ) spzeros(size(rf.M(θ), 1), length(ν)) spdiagm(rf.M(θ) * x - rf.h(θ))
    ]
end