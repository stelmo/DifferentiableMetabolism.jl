
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

Differentiate `m` with respect to `parameters` at solution `sol`, and substitute in `parameter_values`.
"""
function differentiate_model(
    m::C.ConstraintTree,
    sol;
    parameter_values::Dict{Symbol,Float64},
    optimizer,
    zero_tol = 1e-6,
    objective_id = :objective, # must be at top level
    scale = true,
)

    zero_idxs = Set{Int}() # all variables that are close to zero
    C.zip(m, sol) do x, y
        if length(x.value.idxs) == 1 && abs(y) <= zero_tol # TODO should have gene and flux specific zeros?
            push!(zero_idxs, first(x.value.idxs))
        end
        x # TODO this feels a bit wasted?
    end

    vars = [idx in zero_idxs ? zero(C.LinearValue) : C.variable(; idx).value for idx = 1:C.variable_count(m)]

    pm = C.prune_variables(C.substitute(m, vars)) # CTs is wonderful

    pm = C.drop_zeros(pm)

    pm = C.filter_leaves(pm) do leaf
        !isempty(leaf.value.idxs)
    end
    
    parameters = Set{Symbol}()
    C.traverse(pm) do leaf
        ps = []
        leaf.value isa LinearValueP && append!(ps, leaf.value.weights)
        leaf.value isa QuadraticValueP && append!(ps, leaf.value.weights)
        leaf.bound isa EqualToP && push!(ps, leaf.bound.equal_to)
        leaf.bound isa BetweenP && append!(ps, [leaf.bound.lower, leaf.bound.upper])

        for s in filter(!isreal, F.value.(ps))
            push!(parameters, s)
        end
    end

    psol = optimized_constraints_with_parameters(
        pm,
        parameter_values;
        objective = pm[objective_id].value,
        optimizer,
    )

    pm_kkt, vids = differentiate_prepare_kkt(pm, pm[objective_id].value, collect(parameters))

    sens = differentiate_solution(
        pm_kkt,
        psol.primal_values,
        psol.equality_dual_values,
        psol.inequality_dual_values,
        parameter_values,
        scale = scale, 
    )

    sens, vids, collect(parameters)
end

