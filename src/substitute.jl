
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

"""
$(TYPEDSIGNATURES)

Straightforward recursive evaluator for `Ex`s.
"""
substitute(x::Ex, lookup) = substitute(x.node_value, x.children, lookup)
substitute(x, lookup) = substitute(x, nothing, lookup)
substitute(x::Symbol, ::Nothing, lookup) = lookup(x)
substitute(x::Symbol, _, lookup) = lookup(x)
substitute(x, ::Nothing, _) = x
substitute(x, cs, lookup) = x(substitute.(cs, Ref(lookup))...)

#=
Extend the substitute function to convert parameters to numbers.
=#
substitute(x::LinearValueP, lookup) =
    C.LinearValue(x.idxs, substitute.(x.weights, Ref(lookup)))

substitute(x::QuadraticValueP, lookup) =
    C.QuadraticValue(x.idxs, substitute.(x.weights, Ref(lookup)))

substitute(x::BetweenP, lookup) =
    C.Between(substitute(x.lower, lookup), substitute(x.upper, lookup))

substitute(x::EqualToP, lookup) = C.EqualTo(substitute(x.equal_to, lookup))

substitute(x::C.Between, lookup) = x

substitute(x::C.EqualTo, lookup) = x

substitute(x::C.LinearValue, lookup) = x

substitute(x::Nothing, lookup) = nothing

substitute(x::C.Constraint, lookup) =
    C.Constraint(substitute(C.value(x), lookup), substitute(C.bound(x), lookup))

substitute(x::C.ConstraintTree, lookup) = C.map(c -> substitute(c, lookup), x, C.Constraint)
