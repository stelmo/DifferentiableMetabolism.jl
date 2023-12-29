
@kwdef mutable struct ParameterBetween <: C.Bound
    lower::S.Num
    upper::S.Num
end

export ParameterBetween

Base.:-(x::ParameterBetween) = -1 * x
Base.:*(a::ParameterBetween, b::Real) = b * a
Base.:*(a::Real, b::ParameterBetween) = ParameterBetween(a * b.lower, a * b.upper)
Base.:/(a::ParameterBetween, b::Real) = ParameterBetween(a.lower / b, a.upper / b)
