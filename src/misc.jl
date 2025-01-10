
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

C.Constraint(v::T, b::Tuple{X,Y}) where {T<:C.Value,X<:Ex,Y<:Ex} =
    C.Constraint(v, BetweenP(b...))

"""
$(TYPEDSIGNATURES)

Return the names of variables in `m`.

NOTE: this function assumes that all variables are bounded explicitly in the
model. `nothing` bounds are ignored.
"""
function variable_order(m)
    c = []
    ff(p, x::C.ConstraintTree) = nothing
    ff(p, x::C.Constraint) = begin
        if length(x.value.idxs) == 1 && !isnothing(x.bound) #TODO assumes that all variables are bounded somehow!
            idx = first(x.value.idxs)
            push!(c, (idx, p))
        end
    end

    C.itraverse(ff, m)

    _idxs = sortperm(first.(c)) # still need to do this because sets don't respect insertion order
    last.(c)[_idxs]
end
