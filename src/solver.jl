
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

Builds a matrix representation of bounds.
"""
function constraint_matrix_vector(eqs, m, parameters)
    Ib = Int64[]
    Vb = Float64[]
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]

    parameter_lookup(k) = get(parameters, k, error("Parameter $k not supplied."))

    for (i, (val, rhs)) in enumerate(eqs)
        rhs = substitute(rhs, parameter_lookup)
        a = substitute(val, parameter_lookup)

        # TODO this seems prone to error
        if rhs != 0.0 && rhs != -0.0
            push!(Ib, i)
            push!(Vb, rhs)
        end
        append!(Is, fill(i, length(a.idxs)))
        append!(Js, a.idxs)
        append!(Vs, a.weights)
    end
    SA.sparse(Is, Js, Vs, length(eqs), C.var_count(m)), SA.sparsevec(Ib, Vb, length(eqs))
end

"""
$(TYPEDSIGNATURES)

Construct a JuMP model by substituting `parameters` into the model, `m`. Set the
`objective` and the `optimizer`, as well as the `sense` similar to
`COBREXA.optimization_model`.

Converts all inequality constraints to the form `A * x ≤ b`.
"""
function optimization_model_with_parameters(
    m::C.ConstraintTree,
    parameters::Dict{Symbol,Float64};
    objective::C.Value,
    optimizer,
    sense,
)
    model = J.Model(optimizer)

    J.@variable(model, x[1:C.var_count(m)])
    J.@objective(model, sense, C.substitute(substitute(objective, k -> parameters[k]), x))

    eqs = equality_constraints(m)
    ineqs = inequality_constraints(m)
    # variables with nothing bounds are implicitly handled by the solver

    # E * x = d
    E, d = constraint_matrix_vector(eqs, m, parameters)
    J.@constraint(model, eqcons, E * x .== d)

    # M * x ≤ h
    M, h = constraint_matrix_vector(ineqs, m, parameters)
    J.@constraint(model, ineqcons, M * x .<= h)

    return model
end

"""
$(TYPEDSIGNATURES)

Solve a `model` using `optimizer` by substituting in `parameters`. Optional
arguments are the same as in COBREXA.

If the model does not solve successfully return `nothing`. Otherwise, return a
named tuple of the solution tree, and vectors containing the values of the
primal variables, the dual variables.

These duals are ordered according to the constraint output of calling
[`equality_constraints`](@ref) and [`inequality_constraints`](@ref)
respectively.
"""
function optimized_constraints_with_parameters(
    model::C.ConstraintTree,
    parameters::Dict{Symbol,Float64};
    optimizer,
    modifications = [],
    objective::C.Value,
    sense = X.Maximal,
)
    om = optimization_model_with_parameters(model, parameters; objective, optimizer, sense)
    for m in modifications
        m(om)
    end
    J.optimize!(om)

    X.is_solved(om) ?
    (
        tree = C.substitute_values(substitute(m, k -> parameters[k]), J.value.(om[:x])),
        primal_values = J.value.(om[:x]),
        equality_duals = J.dual.(om[:eqcons]),
        inequality_duals = J.dual.(om[:ineqcons]),
    ) : nothing
end

