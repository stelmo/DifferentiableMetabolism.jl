
@kwdef mutable struct ParameterBetween <: C.Bound
    lower::S.Num
    upper::S.Num
end

ParameterBetween(x::Union{Float64,Int}, y::S.Num) = ParameterBetween(S.Num(x), y)
ParameterBetween(x::S.Num, y::Union{Float64,Int}) = ParameterBetween(x, S.Num(y))

export ParameterBetween

Base.:-(x::ParameterBetween) = -1 * x
Base.:*(a::ParameterBetween, b::Real) = b * a
Base.:*(a::Real, b::ParameterBetween) = ParameterBetween(a * b.lower, a * b.upper)
Base.:/(a::ParameterBetween, b::Real) = ParameterBetween(a.lower / b, a.upper / b)


@kwdef mutable struct ParameterEqualTo <: C.Bound
    equal_to::S.Num
end

ParameterEqualTo(y::Union{Float64,Int}) = ParameterBetween(S.Num(y))

export ParameterEqualTo

Base.:-(x::ParameterEqualTo) = -1 * x
Base.:*(a::ParameterEqualTo, b::Real) = b * a
Base.:/(a::ParameterEqualTo, b::Real) = ParameterEqualTo(a.equal_to / b)
Base.:*(a::Real, b::ParameterEqualTo) = ParameterEqualTo(a * b.equal_to)
