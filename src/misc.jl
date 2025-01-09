
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

# this doesn't really fit in anywhere

C.drop_zeros(x::LinearValueP) = LinearValueP(
    idxs = x.idxs[F.value.(x.weights).!=0],
    weights = x.weights[F.value.(x.weights).!=0],
)

C.drop_zeros(x::QuadraticValueP) = QuadraticValueT(
    idxs = x.idxs[F.value.(x.weights).!=0],
    weights = x.weights[F.value.(x.weights).!=0],
)

Base.isreal(x::Symbol) = false

C.Constraint(v::T, b::Ex) where {T<:C.Value} = C.Constraint(v, EqualToP(b))

C.Constraint(v::T, b::Tuple{X,Y}) where {T<:C.Value,X<:Ex,Y<:Ex} = C.Constraint(v, BetweenP(b...))
