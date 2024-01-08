
#=
Copyright (c) 2023, Heinrich-Heine University Duesseldorf
Copyright (c) 2023, University of Luxembourg

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

#=
Extend the substitute function from symbolics to convert parameters to numbers.
Converts all the ParameterXXX types to their XXX types normally found in
ConstraintTrees.
=#

Symbolics.substitute(x::ParameterLinearValue, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.LinearValue(
        x.idxs,
        Symbolics.value.(Symbolics.substitute(x.weights, rule)),
    )

Symbolics.substitute(x::ParameterQuadraticValue, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.QuadraticValue(
        x.idxs,
        Symbolics.value.(Symbolics.substitute(x.weights, rule)),
    )

Symbolics.substitute(x::ParameterBetween, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.Between(
        Symbolics.value.(Symbolics.substitute(x.lower, rule)),
        Symbolics.value.(Symbolics.substitute(x.upper, rule)),
    )

Symbolics.substitute(x::ParameterEqualTo, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.EqualTo(Symbolics.value.(Symbolics.substitute(x.equal_to, rule)))

Symbolics.substitute(x::ConstraintTrees.Between, rule::Dict{Symbolics.Num,Float64}) = x

Symbolics.substitute(x::ConstraintTrees.EqualTo, rule::Dict{Symbolics.Num,Float64}) = x

Symbolics.substitute(x::ConstraintTrees.LinearValue, rule::Dict{Symbolics.Num,Float64}) = x

Symbolics.substitute(x::Nothing, rule::Dict{Symbolics.Num,Float64}) = nothing

Symbolics.substitute(x::ConstraintTrees.Constraint, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.Constraint(
        Symbolics.substitute(ConstraintTrees.value(x), rule),
        Symbolics.substitute(ConstraintTrees.bound(x), rule),
    )

Symbolics.substitute(x::ConstraintTrees.ConstraintTree, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.map(
        c -> Symbolics.substitute(c, rule),
        x,
        ConstraintTrees.Constraint,
    )
