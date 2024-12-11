
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
ConstraintTrees.substitute(x::ParameterLinearValue, y) =
    ConstraintTrees.sum(
        (idx == 0 ? x.weights[i] : x.weights[i] * y[idx] for (i, idx) in enumerate(x.idxs)),
        init = zero(eltype(y)),
    )

ConstraintTrees.substitute(x::ParameterQuadraticValue, y) =
    ConstraintTrees.sum(
        (
            let (idx1, idx2) = x.idxs[i]
                (idx1 == 0 ? 1.0 : y[idx1]) * (idx2 == 0 ? 1.0 : y[idx2]) * w
            end for (i, w) in enumerate(x.weights)
        ),
        init = zero(eltype(y)),
    )

ConstraintTrees.incr_var_idxs(x::ParameterLinearValue, incr::Int) =
    ParameterLinearValue(idxs = ConstraintTrees.incr_var_idx.(x.idxs, incr), weights = x.weights)

ConstraintTrees.incr_var_idxs(x::ParameterQuadraticValue, incr::Int) = ParameterQuadraticValue(
    idxs = broadcast(ii -> ConstraintTrees.incr_var_idx.(ii, incr), x.idxs),
    weights = x.weights,
)

ConstraintTrees.collect_variables!(x::ParameterLinearValue, out) =
    for idx in x.idxs
        push!(out, idx)
    end

ConstraintTrees.collect_variables!(x::ParameterQuadraticValue, out) =
    for (idx, idy) in x.idxs
        push!(out, idx, idy)
    end

ConstraintTrees.renumber_variables(x::ParameterLinearValue, mapping) =
    ParameterLinearValue(idxs = [mapping[idx] for idx in x.idxs], weights = x.weights)

ConstraintTrees.renumber_variables(x::ParameterQuadraticValue, mapping) = 
    ParameterLinearValue(
    idxs = [(mapping[idx], mapping[idy]) for (idx, idy) in x.idxs],
    weights = x.weights,
)
