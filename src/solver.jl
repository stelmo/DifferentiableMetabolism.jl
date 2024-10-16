
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

Changes from copied code are indicated.
=#

# TODO the substitute of symbolics is (unnecessarily) much more powerful than
# the "value with these variable values" that we actually need. Unfortunately
# Symbolics don't realy have anything that would "just do" the simple thing.
# Thus this hack.
fast_subst(x::Symbolics.Num, y) = Symbolics.fast_substitute(x, y)
fast_subst(x, y) = Symbolics.substitute(x, y)

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
    for (i, (val, rhs)) in enumerate(eqs)
        rhs = fast_subst(rhs, parameters)
        a = fast_subst(val, parameters)

        # TODO this seems prone to error
        if Symbolics.value(rhs) != 0.0 && Symbolics.value(rhs) != -0.0
            push!(Ib, i)
            push!(Vb, Symbolics.value(rhs))
        end
        append!(Is, fill(i, length(a.idxs)))
        append!(Js, a.idxs)
        append!(Vs, a.weights)
    end
    SparseArrays.sparse(Is, Js, Vs, length(eqs), ConstraintTrees.var_count(m)),
    SparseArrays.sparsevec(Ib, Vb, length(eqs))
end

"""
$(TYPEDSIGNATURES)

Construct a JuMP model by substituting `parameters` into the model, `m`. Set the
`objective` and the `optimizer`, as well as the `sense` similar to
`COBREXA.optimization_model`.

Converts all inequality constraints to the form `A * x ≤ b`.
"""
function optimization_model_with_parameters(
    m::ConstraintTrees.ConstraintTree,
    parameters::Dict{Symbolics.Num,Float64};
    objective::ConstraintTrees.Value,
    optimizer,
    sense,
)
    model = JuMP.Model(optimizer)

    JuMP.@variable(model, x[1:ConstraintTrees.var_count(m)])
    JuMP.@objective(
        model,
        sense,
        ConstraintTrees.substitute(Symbolics.substitute(objective, parameters), x)
    )

    eqs = equality_constraints(m)
    ineqs = inequality_constraints(m)
    # variables with nothing bounds are implicitly handeled by the solver

    # E * x = d
    E, d = constraint_matrix_vector(eqs, m, parameters)
    JuMP.@constraint(model, eqcons, E * x .== d)

    # M * x ≤ h
    M, h = constraint_matrix_vector(ineqs, m, parameters)
    JuMP.@constraint(model, ineqcons, M * x .<= h)

    return model
end

export optimization_model_with_parameters

"""
$(TYPEDSIGNATURES)

Solve a model, `m`, by forwarding arguments to
[`optimization_model_with_parameters`](@ref). 

Optionally, set optimizer attributes with `modifications`. If the model does not
solve successfully return `nothing`. Otherwise, return a tuple of the solution
tree, and vectors containing the values of the primal variables, the equality
constraint dual variables. 

These duals are ordered according to the constraint output of calling
[`equality_constraints`](@ref) and [`inequality_constraints`](@ref)
respectively.
"""
function optimized_constraints_with_parameters(
    m::ConstraintTrees.ConstraintTree,
    parameters::Dict{Symbolics.Num,Float64};
    modifications = [],
    objective::ConstraintTrees.Value,
    optimizer,
    sense = COBREXA.Maximal,
)
    om = optimization_model_with_parameters(m, parameters; objective, optimizer, sense)
    for m in modifications
        m(om)
    end
    JuMP.optimize!(om)

    COBREXA.is_solved(om) ?
    (
        ConstraintTrees.substitute_values(
            Symbolics.substitute(m, parameters),
            JuMP.value.(om[:x]),
        ),
        JuMP.value.(om[:x]),
        JuMP.dual.(om[:eqcons]),
        JuMP.dual.(om[:ineqcons]),
    ) : nothing
end

export optimized_constraints_with_parameters
