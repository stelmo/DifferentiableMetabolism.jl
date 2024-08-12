
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

#=
Extend various functions in ConstraintTrees to work with their Parameter based equivalents.
=#

ConstraintTrees.var_count(x::ParameterLinearValue) = isempty(x.idxs) ? 0 : last(x.idxs)

ConstraintTrees.var_count(x::ParameterQuadraticValue) =
    isempty(x.idxs) ? 0 : let (_, max) = last(x.idxs)
        max
    end

# substitute in the variables as symbolic numbers - useful to construct the KKT function
ConstraintTrees.substitute(x::ParameterLinearValue, y::Vector{Symbolics.Num}) = sum(
    (idx == 0 ? x.weights[i] : x.weights[i] * y[idx] for (i, idx) in enumerate(x.idxs)),
    init = Symbolics.Num(0.0),
)

ConstraintTrees.substitute(x::ParameterQuadraticValue, y::Vector{Symbolics.Num}) = sum(
    (
        let (idx1, idx2) = x.idxs[i]
            (idx1 == 0 ? 1.0 : y[idx1]) * (idx2 == 0 ? 1.0 : y[idx2]) * w
        end for (i, w) in enumerate(x.weights)
    ),
    init = Symbolics.Num(0.0),
)
