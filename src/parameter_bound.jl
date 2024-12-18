
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
=#

"""
$(TYPEDEF)

Representation of an "interval" bound where the lower and upper bound values are
parameters. Since `Expression` is a subtype of `Real`, the bounds could also
be any real number, but they are converted by the constructors to
`Expression`s.

# Fields
$(TYPEDFIELDS)
"""
@kwdef mutable struct ParameterBetween <: ConstraintTrees.Bound
    lower::Expression
    upper::Expression
end

# ParameterBetween(x::Union{Float64,Int,Expression}, y::Union{Float64,Int,Expression}) =
#     ParameterBetween(Expression(x), Expression(y))

export ParameterBetween

Base.:-(x::ParameterBetween) = -1 * x
Base.:*(a::ParameterBetween, b::Real) = b * a
Base.:/(a::ParameterBetween, b::Real) = (1 / b) * a
Base.:*(a::Real, b::ParameterBetween) = ParameterBetween(a * b.lower, a * b.upper)

"""
$(TYPEDEF)

Representation of an "equality" bound, where the bound value is a parameter.
Since `Expression` is a subtype of `Real`, the bound could also be any real
number, but it is converted by the constructor to a `Expression`.

# Fields
$(TYPEDFIELDS)
"""
@kwdef mutable struct ParameterEqualTo <: ConstraintTrees.Bound
    equal_to::Expression
end

ParameterEqualTo(y::Union{Float64,Int}) = ParameterEqualTo(Expression(y))

export ParameterEqualTo

Base.:-(x::ParameterEqualTo) = -1 * x
Base.:*(a::ParameterEqualTo, b::Real) = b * a
Base.:/(a::ParameterEqualTo, b::Real) = (1 / b) * a
Base.:*(a::Real, b::ParameterEqualTo) = ParameterEqualTo(a * b.equal_to)
