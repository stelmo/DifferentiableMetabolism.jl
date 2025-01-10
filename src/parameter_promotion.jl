
#=
Copyright (c) 2025, Heinrich-Heine University Duesseldorf
Copyright (c) 2025, University of Luxembourg

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


C.LinearValueT(x::Ex) =
    F.is_identically_zero(x) ? LinearValueP(; idxs = Int[], weights = Ex[]) :
    LinearValueP(; idxs = [0], weights = [x])
C.QuadraticValueT(x::Ex) =
    F.is_identically_zero(x) ?
    QuadraticValueP(; idxs = Vector{Tuple{Int,Int}}(), weights = Ex[]) :
    QuadraticValueP(; idxs = [(0, 0)], weights = [x])

## Promotion

# Exs and LinearValues
Base.:+(a::Ex, b::C.LinearValue) = b + a
Base.:+(a::C.LinearValue, b::Ex) = a + C.LinearValueT(b)

Base.:-(a::Ex, b::C.LinearValue) = -b + a
Base.:-(a::C.LinearValue, b::Ex) = a - C.LinearValueT(b)

Base.:*(a::Ex, b::C.LinearValue) = b * a
Base.:*(a::C.LinearValue, b::Ex) = C.LinearValueT(a.idxs, b .* a.weights)

Base.:/(a::C.LinearValue, b::Ex) = C.LinearValueT(a.idxs, a.weights ./ b)

# LinearValue and LinearValuePs
Base.:+(a::C.LinearValue, b::LinearValueP) = b + a
Base.:+(a::LinearValueP, b::C.LinearValue) = a + C.LinearValueT(b.idxs, Ex.(b.weights))

Base.:-(a::C.LinearValue, b::LinearValueP) = -b + a
Base.:-(a::LinearValueP, b::C.LinearValue) = a - C.LinearValueT(b.idxs, Ex.(b.weights))

# # Exs and QuadraticValueTs
Base.:+(a::Ex, b::C.QuadraticValue) = b + a
Base.:+(a::C.QuadraticValue, b::Ex) = a + C.QuadraticValueT(b)

Base.:-(a::Ex, b::C.QuadraticValue) = -b + a
Base.:-(a::C.QuadraticValue, b::Ex) = a - C.QuadraticValueT(b)

Base.:*(a::Ex, b::C.QuadraticValue) = b * a
Base.:*(a::C.QuadraticValue, b::Ex) = C.QuadraticValueT(a.idxs, b .* a.weights)

Base.:/(a::C.QuadraticValue, b::Ex) = C.QuadraticValueT(a.idxs, a.weights ./ b)

# QuadraticValues and QuadraticValueTs
Base.:+(a::C.QuadraticValue, b::QuadraticValueP) = b + a
Base.:+(a::QuadraticValueP, b::C.QuadraticValue) =
    a + C.QuadraticValueT(b.idxs, Ex.(b.weights))

Base.:-(a::C.QuadraticValue, b::QuadraticValueP) = -b + a
Base.:-(a::QuadraticValueP, b::C.QuadraticValue) =
    a - C.QuadraticValueT(b.idxs, Ex.(b.weights))

## Cross interaction terms

# LinearValues and QuadraticValueTs (+, - only)
Base.:+(a::C.LinearValue, b::QuadraticValueP) = b + a
Base.:+(a::QuadraticValueP, b::C.LinearValue) = a + C.LinearValueT(b.idxs, Ex.(b.weights))

Base.:-(a::C.LinearValue, b::QuadraticValueP) = -b + a
Base.:-(a::QuadraticValueP, b::C.LinearValue) = a - C.LinearValueT(b.idxs, Ex.(b.weights))

# LinearValue * LinearValueP
Base.:*(a::LinearValueP, b::C.LinearValue) = b * a
Base.:*(a::C.LinearValue, b::LinearValueP) = C.LinearValueT(a.idxs, Ex.(a.weights)) * b

# downcast
LinearValueP(x::Float64) = C.LinearValue(x)
LinearValueP(x::Int64) = C.LinearValue(x)
QuadraticValueP(x::Float64) = C.QuadraticValue(x)
QuadraticValueP(x::Int64) = C.QuadraticValue(x)

# Expressions and bounds
Base.:*(a::Ex, b::C.Between) = b * a
Base.:*(a::C.Between, b::Ex) = BetweenP(; lower = a.lower * b, upper = a.upper * b)
Base.:/(a::C.Between, b::Ex) = BetweenP(; lower = a.lower / b, upper = a.upper / b)

Base.:*(a::Ex, b::C.EqualTo) = b * a
Base.:*(a::C.EqualTo, b::Ex) = EqualToP(; equal_to = a.equal_to * b)
Base.:/(a::C.EqualTo, b::Ex) = EqualToP(; equal_to = a.equal_to / b)
