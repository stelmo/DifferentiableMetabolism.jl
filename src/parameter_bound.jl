
"""
$(TYPEDEF)

Representation of an "interval" bound where the lower and upper bound values are
parameters. Since [`Symbolics.Num`](@ref) is a subtype of `Real`, the bounds
could also be any real number, but they are converted by the constructors to
[`Symbolics.Num`](@ref)s. 

# Fields
$(TYPEDFIELDS)
"""
@kwdef mutable struct ParameterBetween <: ConstraintTrees.Bound
    lower::Symbolics.Num
    upper::Symbolics.Num
end

ParameterBetween(x::Union{Float64,Int,Symbolics.Num}, y::Union{Float64,Int,Symbolics.Num}) =
    ParameterBetween(Symbolics.Num(x), Symbolics.Num(y))

export ParameterBetween

Base.:-(x::ParameterBetween) = -1 * x
Base.:*(a::ParameterBetween, b::Real) = b * a
Base.:*(a::Real, b::ParameterBetween) = ParameterBetween(a * b.lower, a * b.upper)
Base.:/(a::ParameterBetween, b::Real) = ParameterBetween(a.lower / b, a.upper / b)

"""
$(TYPEDEF)

Representation of an "equality" bound, where the bound value is a parameter.
Since [`Symbolics.Num`](@ref) is a subtype of `Real`, the bound could also be
any real number, but it is converted by the constructor to a
[`Symbolics.Num`](@ref).

# Fields
$(TYPEDFIELDS)
"""
@kwdef mutable struct ParameterEqualTo <: ConstraintTrees.Bound
    equal_to::Symbolics.Num
end

ParameterEqualTo(y::Union{Float64,Int}) = ParameterBetween(Symbolics.Num(y))

export ParameterEqualTo

Base.:-(x::ParameterEqualTo) = -1 * x
Base.:*(a::ParameterEqualTo, b::Real) = b * a
Base.:/(a::ParameterEqualTo, b::Real) = ParameterEqualTo(a.equal_to / b)
Base.:*(a::Real, b::ParameterEqualTo) = ParameterEqualTo(a * b.equal_to)
