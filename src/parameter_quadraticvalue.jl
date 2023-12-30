
"""
$(TYPEDEF)

An extension of [`ConstraintTrees.QuadraticValue`](@ref) where the weights are
parameters.

Behaves similarly to [`ConstraintTrees.QuadraticValue`](@ref). Thus, the
cleanest way to construct a `ParameterQuadraticValue` is to multiply
two [`ParameterLinearValue`](@ref)s.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef struct ParameterQuadraticValue <: ConstraintTrees.Value
    idxs::Vector{Tuple{Int,Int}}
    weights::Vector{Symbolics.Num}
end

"""
$(TYPEDSIGNATURES)

Construct a constant [`ParameterQuadraticValue`](@ref) with a single affine element.
"""
ParameterQuadraticValue(x::Real) =
    iszero(x) ? ParameterQuadraticValue(idxs = [], weights = []) :
    ParameterQuadraticValue(idxs = [(0, 0)], weights = [Symbolics.Num(x)])

"""
$(TYPEDSIGNATURES)

Construct a [`ParameterQuadraticValue`](@ref) that is equivalent to a given [`LinearValue`](@ref).
"""
ParameterQuadraticValue(x::LinearValue) =
    ParameterQuadraticValue(idxs = [(0, idx) for idx in x.idxs], weights = x.weights)

Base.convert(::Type{ParameterQuadraticValue}, x::Real) = ParameterQuadraticValue(x)
Base.convert(::Type{ParameterQuadraticValue}, x::LinearValue) = ParameterQuadraticValue(x)
Base.zero(::Type{ParameterQuadraticValue}) =
    ParameterQuadraticValue(idxs = [], weights = [])
Base.:+(a::Real, b::ParameterQuadraticValue) = ParameterQuadraticValue(a) + b
Base.:+(a::ParameterQuadraticValue, b::Real) = a + ParameterQuadraticValue(b)
Base.:+(a::LinearValue, b::ParameterQuadraticValue) = ParameterQuadraticValue(a) + b
Base.:+(a::ParameterQuadraticValue, b::LinearValue) = a + ParameterQuadraticValue(b)
Base.:-(a::ParameterQuadraticValue) = -1 * a
Base.:-(a::Real, b::ParameterQuadraticValue) = ParameterQuadraticValue(a) - b
Base.:-(a::ParameterQuadraticValue, b::Real) = a - ParameterQuadraticValue(b)
Base.:-(a::LinearValue, b::ParameterQuadraticValue) = ParameterQuadraticValue(a) - b
Base.:-(a::ParameterQuadraticValue, b::LinearValue) = a - ParameterQuadraticValue(b)
Base.:*(a::Real, b::ParameterQuadraticValue) = b * a
Base.:*(a::ParameterQuadraticValue, b::Real) =
    ParameterQuadraticValue(idxs = a.idxs, weights = b .* a.weights)
Base.:-(a::ParameterQuadraticValue, b::ParameterQuadraticValue) = a + (-1 * b)
Base.:/(a::ParameterQuadraticValue, b::Real) =
    ParameterQuadraticValue(idxs = a.idxs, weights = a.weights ./ b)

"""
$(TYPEDSIGNATURES)

Internal helper for co-lex ordering of indexes.
"""
colex_le((a, b), (c, d)) = (b, a) < (d, c)

function Base.:+(a::ParameterQuadraticValue, b::ParameterQuadraticValue)
    r_idxs = Tuple{Int,Int}[]
    r_weights = Float64[]
    ai = 1
    ae = length(a.idxs)
    bi = 1
    be = length(b.idxs)

    while ai <= ae && bi <= be
        if colex_le(a.idxs[ai], b.idxs[bi])
            push!(r_idxs, a.idxs[ai])
            push!(r_weights, a.weights[ai])
            ai += 1
        elseif colex_le(b.idxs[bi], a.idxs[ai])
            push!(r_idxs, b.idxs[bi])
            push!(r_weights, b.weights[bi])
            bi += 1
        else # index pairs are equal; merge case
            push!(r_idxs, a.idxs[ai])
            push!(r_weights, a.weights[ai] + b.weights[bi])
            ai += 1
            bi += 1
        end
    end
    while ai <= ae
        push!(r_idxs, a.idxs[ai])
        push!(r_weights, a.weights[ai])
        ai += 1
    end
    while bi <= be
        push!(r_idxs, b.idxs[bi])
        push!(r_weights, b.weights[bi])
        bi += 1
    end
    ParameterQuadraticValue(idxs = r_idxs, weights = r_weights)
end

Base.:*(a::LinearValue, b::LinearValue) =
    let vals = a.weights .* b.weights'
        ParameterQuadraticValue(
            idxs = [(aidx, bidx) for bidx in b.idxs for aidx in a.idxs if aidx <= bidx],
            weights = [
                vals[ai, bi] for bi in eachindex(b.idxs) for
                ai in eachindex(a.idxs) if a.idxs[ai] <= b.idxs[bi]
            ],
        ) + ParameterQuadraticValue(
            idxs = [(bidx, aidx) for aidx in a.idxs for bidx in b.idxs if bidx < aidx],
            weights = [
                vals[ai, bi] for ai in eachindex(a.idxs) for
                bi in eachindex(b.idxs) if b.idxs[bi] < a.idxs[ai]
            ],
        )
    end

"""
$(TYPEDSIGNATURES)

Substitute anything vector-like as variable values into the [`ParameterQuadraticValue`](@ref)
and return the result.
"""
substitute(x::ParameterQuadraticValue, y) = sum(
    (
        let (idx1, idx2) = x.idxs[i]
            (idx1 == 0 ? 1.0 : y[idx1]) * (idx2 == 0 ? 1.0 : y[idx2]) * w
        end for (i, w) in enumerate(x.weights)
    ),
    init = 0.0,
)


