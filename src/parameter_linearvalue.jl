
"""
$(TYPEDEF)

An extension of [`ConstraintTrees.LinearValue`](@ref) where the weights are parameters.

`ParameterLinearValue`s can be combined additively and multiplied by real-number constants.

Multiplying two `ParameterLinearValue`s yields a quadratic form (in a [`ParameterQuadraticValue`](@ref)).

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct ParameterLinearValue <: ConstraintTrees.Value
    idxs::Vector{Int}
    weights::Vector{Symbolics.Num}
end

# convert all weights to Nums
ParameterLinearValue(idxs::Vector{Int}, weights::Vector{Union{Int,Float64}}) =
    ParameterLinearValue(idxs, convert.(Symbolics.Num, weights))

export ParameterLinearValue

ParameterLinearValue(x::Real) =
    iszero(x) ? ParameterLinearValue(idxs = [], weights = []) :
    ParameterLinearValue(idxs = [0], weights = [Symbolics.Num(x)])

Base.convert(::Type{ParameterLinearValue}, x::Real) = ParameterLinearValue(x)
Base.zero(::Type{ParameterLinearValue}) = ParameterLinearValue(idxs = [], weights = [])
Base.:+(a::Real, b::ParameterLinearValue) = ParameterLinearValue(a) + b
Base.:+(a::ParameterLinearValue, b::Real) = a + ParameterLinearValue(b)
Base.:-(a::Real, b::ParameterLinearValue) = ParameterLinearValue(a) - b
Base.:-(a::ParameterLinearValue, b::Real) = a - ParameterLinearValue(b)
Base.:-(a::ParameterLinearValue, b::ParameterLinearValue) = a + (-1 * b)
Base.:-(a::ParameterLinearValue) = -1 * a
Base.:*(a::Real, b::ParameterLinearValue) = b * a
Base.:*(a::ParameterLinearValue, b::Real) = ParameterLinearValue(a.idxs, b .* a.weights)
Base.:/(a::ParameterLinearValue, b::Real) = ParameterLinearValue(a.idxs, a.weights ./ b)

# Promote LinearValue to ParameterLinearValue
Base.:+(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = ParameterLinearValue(a) + b
Base.:+(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = a + ParameterLinearValue(b)
Base.:-(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = ParameterLinearValue(a) - b
Base.:-(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = a - ParameterLinearValue(b)
Base.:*(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = b * a
Base.:*(a::ConstraintTrees.LinearValue, b::Symbolics.Num) =
    ParameterLinearValue(a.idxs, b .* a.weights)
Base.:/(a::ConstraintTrees.LinearValue, b::Symbolics.Num) =
    ParameterLinearValue(a.idxs, a.weights ./ b)

Base.:+(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) =
    ParameterLinearValue(a.idxs, a.weights) + b
Base.:+(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) =
    a + ParameterLinearValue(b.idxs, b.weights)
Base.:-(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) =
    ParameterLinearValue(a.idxs, a.weights) - b
Base.:-(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) =
    a - ParameterLinearValue(b.idxs, b.weights)

# add two ParameterLinearValues
function Base.:+(a::ParameterLinearValue, b::ParameterLinearValue)
    r_idxs = Int[]
    r_weights = Symbolics.Num[] # supertype including Num and Float64
    ai = 1
    ae = length(a.idxs)
    bi = 1
    be = length(b.idxs)
    while ai <= ae && bi <= be
        if a.idxs[ai] < b.idxs[bi]
            push!(r_idxs, a.idxs[ai])
            push!(r_weights, a.weights[ai])
            ai += 1
        elseif a.idxs[ai] > b.idxs[bi]
            push!(r_idxs, b.idxs[bi])
            push!(r_weights, b.weights[bi])
            bi += 1
        else # a.idxs[ai] == b.idxs[bi] -- merge case
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
    ParameterLinearValue(idxs = r_idxs, weights = r_weights)
end
