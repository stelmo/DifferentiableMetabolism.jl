
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

const FloatInt = Union{Float64, Int} # TODO ugly
const LinearValueP = LinearValueT{Expression}
const QuadraticValueP = QuadraticValueT{Expression}

LinearValueT(x::Expression) = FastDifferentiation.is_identically_zero(x) ? LinearValueT{Expression}(;idxs=Int[], weights=Expression[]) : LinearValueT(;idxs=[0,], weights=[x])
QuadraticValueT(x::Expression) = FastDifferentiation.is_identically_zero(x) ? QuadraticValueT{Expression}(;idxs=Vector{Tuple{Int,Int}}(), weights=Expression[]) : QuadraticValueT(;idxs=[(0,0),], weights=[x])

# Expressions and LinearValues
Base.:+(a::Expression, b::LinearValue) = b + a
Base.:+(a::LinearValueT{T}, b::Expression) where {T<:FloatInt} = a + LinearValueT(b)

Base.:-(a::Expression, b::LinearValue) = -b + a
Base.:-(a::LinearValueT{T}, b::Expression) where {T<:FloatInt} = a - LinearValueT(b)

Base.:*(a::Expression, b::LinearValue) = b * a
Base.:*(a::LinearValueT{T}, b::Expression) where {T<:FloatInt} = LinearValueT(a.idxs, b .* a.weights)

Base.:/(a::LinearValueT{T}, b::Expression) where {T<:FloatInt} = LinearValueT(a.idxs, a.weights ./ b)

# LinearValue and LinearValueT{Expression}s
Base.:+(a::LinearValueT{T}, b::LinearValueT{Expression}) where {T<:FloatInt} = b + a
Base.:+(a::LinearValueT{Expression}, b::LinearValueT{T}) where {T<:FloatInt} = a + LinearValueT(b.idxs, Expression.(b.weights))

Base.:-(a::LinearValueT{T}, b::LinearValueT{Expression}) where {T<:FloatInt} = -b + a
Base.:-(a::LinearValueT{Expression}, b::LinearValueT{T}) where {T<:FloatInt} = a - LinearValueT(b.idxs, Expression(b.weights))

# # Expressions and QuadraticValueTs
Base.:+(a::Expression, b::QuadraticValueT{T}) where {T<:FloatInt} = b + a
Base.:+(a::QuadraticValueT{T}, b::Expression) where {T<:FloatInt} = a + QuadraticValueT(b)

Base.:-(a::Expression, b::QuadraticValueT{T}) where {T<:FloatInt} = -b + a
Base.:-(a::QuadraticValueT{T}, b::Expression) where {T<:FloatInt} = a - QuadraticValueT(b)

Base.:*(a::Expression, b::QuadraticValueT{T}) where {T<:FloatInt} = b * a
Base.:*(a::QuadraticValueT{T}, b::Expression) where {T<:FloatInt} = QuadraticValueT(a.idxs, b .* a.weights)

Base.:/(a::QuadraticValueT{T}, b::Expression) where {T<:FloatInt} = QuadraticValueT(a.idxs, a.weights ./ b)

# QuadraticValues and QuadraticValueTs
Base.:+(a::QuadraticValueT{T}, b::QuadraticValueT{Expression}) where {T<:FloatInt} = b + a
Base.:+(a::QuadraticValueT{Expression}, b::QuadraticValueT{T}) where {T<:FloatInt} = a + QuadraticValueT(b.idxs, Expression.(b.weights))

Base.:-(a::QuadraticValueT{T}, b::QuadraticValueT{Expression}) where {T<:FloatInt} = -b + a
Base.:-(a::QuadraticValueT{Expression}, b::QuadraticValueT{T}) where {T<:FloatInt} = a - QuadraticValueT(b.idxs, Expression.(b.weights))

# Cross interaction terms

# LinearValues and QuadraticValueTs (+, - only)
Base.:+(a::LinearValueT{T}, b::QuadraticValueT{Expression}) where {T<:FloatInt} = b + a
Base.:+(a::QuadraticValueT{Expression}, b::LinearValueT{T}) where {T<:FloatInt} = a + LinearValueT(b.idxs, Expression.(b.weights))

Base.:-(a::LinearValueT{T}, b::QuadraticValueT{Expression}) where {T<:FloatInt} = -b + a
Base.:-(a::QuadraticValueT{Expression}, b::LinearValueT{T}) where {T<:FloatInt} = a - LinearValueT(b.idxs, Expression.(b.weights))

# LinearValue * LinearValueT{Expression}
Base.:*(a::LinearValueT{Expression}, b::LinearValueT{T}) where {T<:FloatInt} = b * a
Base.:*(a::LinearValueT{T}, b::LinearValueT{Expression}) where {T<:FloatInt} = LinearValueT{Expression}(a.idxs, Expression.(a.weights)) * b
