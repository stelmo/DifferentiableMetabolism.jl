
#=
Copyright (c) 2024, Heinrich-Heine University Duesseldorf
Copyright (c) 2024, University of Luxembourg

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

An extension of `ConstraintTrees.LinearValue` where the weights are parameters.

`ParameterLinearValue`s can be combined additively and multiplied by real-number
constants.

Multiplying two `ParameterLinearValue`s yields a quadratic form (in a
[`ParameterQuadraticValue`](@ref)).

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct ParameterLinearValue <: ConstraintTrees.Value
    idxs::Vector{Int}
    weights::Vector{Expression}
end

# convert all weights to Nums
ParameterLinearValue(idxs::Vector{Int}, weights::Vector{Union{Int,Float64}}) =
    ParameterLinearValue(idxs, convert.(Expression, weights))

ParameterLinearValue(x::ConstraintTrees.LinearValue) =
    ParameterLinearValue(x.idxs, convert.(Expression, x.weights))

export ParameterLinearValue

ParameterLinearValue(x::Real) =
    iszero(x) ? ParameterLinearValue(idxs = [], weights = []) :
    ParameterLinearValue(idxs = [0], weights = [Expression(x)])

Base.convert(::Type{ParameterLinearValue}, x::Real) = ParameterLinearValue(x)

Base.zero(::Type{ParameterLinearValue}) = ParameterLinearValue(idxs = [], weights = [])

Base.:+(a::Real, b::ParameterLinearValue) = b + a

Base.:+(a::ParameterLinearValue, b::Real) = a + ParameterLinearValue(b)

Base.:-(a::Real, b::ParameterLinearValue) = -b + a

Base.:-(a::ParameterLinearValue, b::Real) = a - ParameterLinearValue(b)

Base.:-(a::ParameterLinearValue, b::ParameterLinearValue) = a + (-1 * b)

Base.:-(a::ParameterLinearValue) = -1 * a

Base.:*(a::Real, b::ParameterLinearValue) = b * a

Base.:*(a::ParameterLinearValue, b::Real) = ParameterLinearValue(a.idxs, b .* a.weights)

Base.:/(a::ParameterLinearValue, b::Real) = (1 / b) * a

# add two ParameterLinearValues
function Base.:+(a::ParameterLinearValue, b::ParameterLinearValue)
    # Code mostly copied from ConstraintTrees.jl, but marginally changed some
    # types
    r_idxs = Int[]
    r_weights = Expression[]
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
