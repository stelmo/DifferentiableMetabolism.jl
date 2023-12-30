
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

"""
$(TYPEDSIGNATURES)

Construct a JuMP model by substituting `parameters` into the model, `m`. Set the
`objective` and the `optimizer`, as well as the `sense` similar to
[`COBREXA.optimization_model`](@ref).

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

    # E * x = d
    Ib = Int64[]
    Vb = Float64[]
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for (i, (val, rhs)) in enumerate(eqs)
        rhs = Symbolics.substitute(rhs, parameters)
        a = Symbolics.substitute(val, parameters)

        if Symbolics.value(rhs) != 0
            push!(Ib, i)
            push!(Vb, Symbolics.value(rhs))
        end
        append!(Is, fill(i, length(a.idxs)))
        append!(Js, a.idxs)
        append!(Vs, a.weights)
    end
    A = sparse(Is, Js, Vs, length(eqs), ConstraintTrees.var_count(m))
    b = sparsevec(Ib, Vb, length(eqs))
    JuMP.@constraint(model, eqcons, A * x .== b)

    # M * x ≤ h
    Ih = Int64[]
    Vh = Float64[]
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    k = 0
    for (val, lower, upper) in ineqs
        lower = Symbolics.substitute(lower, parameters)
        upper = Symbolics.substitute(upper, parameters)
        a = Symbolics.substitute(val, parameters)

        # lower: l ≤ x => -x ≤ -l
        if !isinf(Symbolics.value(lower))
            k += 1
            if Symbolics.value(lower) != 0
                push!(Ih, k)
                push!(Vh, -Symbolics.value(lower))
            end
            append!(Is, fill(k, length(a.idxs)))
            append!(Js, a.idxs)
            append!(Vs, -a.weights)
        end

        # upper: x ≤ u
        if !isinf(Symbolics.value(upper))
            k += 1
            if Symbolics.value(upper) != 0
                push!(Ih, k)
                push!(Vh, Symbolics.value(upper))
            end
            append!(Is, fill(k, length(a.idxs)))
            append!(Js, a.idxs)
            append!(Vs, a.weights)
        end
    end
    M = sparse(Is, Js, Vs, k, ConstraintTrees.var_count(m))
    h = sparsevec(Ih, Vh, k)
    JuMP.@constraint(model, ineqcons, M * x .<= h)

    return model
end

export optimization_model_with_parameters

"""
$(TYPEDSIGNATURES)

Solve a model, `m`, by forwarding arguments to
[`optimization_model_with_parameters`]. 

Optionally, set optimizer attributes with `modifications`. If the model does not
solve successfully return `nothing`. Otherwise, if `primal_duals_only` is true,
return a tuple of vectors containing the values of the primal variables, the
equality constraint dual variables, and finally the inequality constraint dual
variables. These duals are ordered according to the constraint output of calling
[`equality_constraints`](@ref) and [`inequality_constraints`](@ref)
respectively. If `primal_duals_only` is false, return a constraint tree with the
solution primal values.
"""
function optimized_constraints_with_parameters(
    m::ConstraintTrees.ConstraintTree,
    parameters::Dict{Symbolics.Num,Float64};
    modifications = [],
    objective::ConstraintTrees.Value,
    optimizer,
    sense = COBREXA.Maximal,
    primal_duals_only = false,
)
    om = optimization_model_with_parameters(m, parameters; objective, optimizer, sense)
    for m in modifications
        m(om)
    end
    JuMP.optimize!(om)

    COBREXA.is_solved(om) ?
    (
        primal_duals_only ?
        (JuMP.value.(om[:x]), JuMP.dual.(om[:eqcons]), JuMP.dual.(om[:ineqcons])) :
        ConstraintTrees.constraint_values(
            Symbolics.substitute(m, parameters),
            JuMP.value.(om[:x]),
        )
    ) : nothing
end

export optimized_constraints_with_parameters
