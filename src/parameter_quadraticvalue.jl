
#=
Copyright (c) 2023, Heinrich-Heine University Duesseldorf
Copyright (c) 2023, University of Luxembourg

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Changes from copied code are indicated.
=#

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

ParameterQuadraticValue(idxs::Vector{Tuple{Int,Int}}, weights::Vector{Union{Int,Float64}}) =
    ParameterQuadraticValue(idxs, convert.(Symbolics.Num, weights))

ParameterQuadraticValue(x::ConstraintTrees.QuadraticValue) =
    ParameterQuadraticValue(x.idxs, convert.(Symbolics.Num, x.weights))

ParameterQuadraticValue(x::ParameterLinearValue) =
    ParameterQuadraticValue(idxs = [(0, idx) for idx in x.idxs], weights = x.weights)

ParameterQuadraticValue(x::Real) =
    iszero(x) ? ParameterQuadraticValue(idxs = [], weights = []) :
    ParameterQuadraticValue(idxs = [(0, 0)], weights = [Symbolics.Num(x)])

Base.convert(::Type{ParameterQuadraticValue}, x::Real) = ParameterQuadraticValue(x)

Base.convert(::Type{ParameterQuadraticValue}, x::ParameterLinearValue) =
    ParameterQuadraticValue(x)

Base.zero(::Type{ParameterQuadraticValue}) =
    ParameterQuadraticValue(idxs = [], weights = [])

Base.:+(a::Real, b::ParameterQuadraticValue) = ParameterQuadraticValue(a) + b

Base.:+(a::ParameterQuadraticValue, b::Real) = a + ParameterQuadraticValue(b)

Base.:+(a::ParameterLinearValue, b::ParameterQuadraticValue) =
    ParameterQuadraticValue(a) + b

Base.:+(a::ParameterQuadraticValue, b::ParameterLinearValue) =
    a + ParameterQuadraticValue(b)

Base.:-(a::ParameterQuadraticValue) = -1 * a

Base.:-(a::Real, b::ParameterQuadraticValue) = ParameterQuadraticValue(a) - b

Base.:-(a::ParameterQuadraticValue, b::Real) = a - ParameterQuadraticValue(b)

Base.:-(a::ParameterLinearValue, b::ParameterQuadraticValue) =
    ParameterQuadraticValue(a) - b

Base.:-(a::ParameterQuadraticValue, b::ParameterLinearValue) =
    a - ParameterQuadraticValue(b)

Base.:*(a::Real, b::ParameterQuadraticValue) = b * a

Base.:*(a::ParameterQuadraticValue, b::Real) =
    ParameterQuadraticValue(idxs = a.idxs, weights = b .* a.weights)

Base.:-(a::ParameterQuadraticValue, b::ParameterQuadraticValue) = a + (-1 * b)

Base.:/(a::ParameterQuadraticValue, b::Real) =
    ParameterQuadraticValue(idxs = a.idxs, weights = a.weights ./ b)

function Base.:+(a::ParameterQuadraticValue, b::ParameterQuadraticValue)
    # Code mostly copied from ConstraintTrees.jl, but marginally changed some
    # types
    r_idxs = Tuple{Int,Int}[]
    r_weights = Symbolics.Num[]
    ai = 1
    ae = length(a.idxs)
    bi = 1
    be = length(b.idxs)

    while ai <= ae && bi <= be
        if ConstraintTrees.colex_le(a.idxs[ai], b.idxs[bi])
            push!(r_idxs, a.idxs[ai])
            push!(r_weights, a.weights[ai])
            ai += 1
        elseif ConstraintTrees.colex_le(b.idxs[bi], a.idxs[ai])
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

Base.:*(a::ParameterLinearValue, b::ParameterLinearValue) =
# Code mostly copied from ConstraintTrees.jl, but marginally changed some
# types
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
