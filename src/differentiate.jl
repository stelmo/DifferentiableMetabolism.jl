
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

Return all the equality constraints of `m` as a tuple `({Parameter}LinearValue,
value)` representing `{P}LV == value` for each entry.
"""
function differentiate(
    m::ConstraintTrees.ConstraintTree,
    objective::ConstraintTrees.Value,
    x_vals::Vector{Float64},
    eq_dual_vals::Vector{Float64},
    ineq_dual_vals::Vector{Float64},
    parameter_values::Dict{Symbolics.Num, Float64},
    parameters::Vector{Symbolics.Num}; # might not diff wrt all params
    zero_tol = 1e-8,
)

    Symbolics.@variables x[1:ConstraintTrees.var_count(m)] # primal
    xs = collect(x)

    nonzero_primal_idxs = [i for i in eachindex(xs) if abs(x_vals[i]) > zero_tol]
    zero_primal_idxs = setdiff(1:ConstraintTrees.var_count(m), nonzero_primal_idxs)
    xzeros = Dict(xs[zero_primal_idxs] .=> 0.0)

    # objective
    f = ConstraintTrees.substitute(objective, xs)
    f = Symbolics.substitute(f, xzeros)

    # equality constraints
    # E * x - b = H = 0
    iH = [
        (i, Symbolics.substitute(ConstraintTrees.substitute(lhs, xs) - rhs, xzeros)) for (i, (lhs, rhs)) in enumerate(equality_constraints(m))
    ]
    # filter out all equations that only contain the zero primals
    filter!(x -> !isequal(last(x), -0.0) && !isequal(last(x), 0.0), iH)
    H = last.(iH)
    eq_dual_idxs = first.(iH)
    
    # inequality constraints (must be built the same as in solver.jl)
    # M * x - h = G ≤ 0
    iG = Vector{Tuple{Int, Symbolics.Num}}()
    i = 0
    for (lhs, lower, upper) in inequality_constraints(m)

        # lower: l ≤ x => -x ≤ -l
        l = Symbolics.value(lower)
        if l isa Float64 && isinf(l)
            nothing
        else
            i += 1
            push!(iG, (i, Symbolics.substitute(-ConstraintTrees.substitute(lhs, xs) + lower, xzeros)))
        end

        # upper: x ≤ u
        u = Symbolics.value(upper)
        if u isa Float64 && isinf(u)
            nothing
        else
            i += 1
            push!(iG, (i, Symbolics.substitute(ConstraintTrees.substitute(lhs, xs) - upper, xzeros)))
        end
    end

    # filter out all equations that only contain the zero primals
    filter!(x -> !isequal(last(x), -0.0) && !isequal(last(x), 0.0), iG)
    G = last.(iG)
    ineq_dual_idxs = first.(iG)

    Symbolics.@variables eq_duals[1:length(H)] ineq_duals[1:length(G)]

    kkt_eqns = [ # negatives because of KKT formulation in JuMP
        Symbolics.sparsejacobian([f], xs[nonzero_primal_idxs])' - Symbolics.sparsejacobian(H, xs[nonzero_primal_idxs])' * eq_duals -
        Symbolics.sparsejacobian(G, xs[nonzero_primal_idxs])' * ineq_duals
        H
        G .* ineq_duals
    ]

    A = Symbolics.sparsejacobian(
        kkt_eqns[:],
        [xs[nonzero_primal_idxs]; eq_duals; ineq_duals],
    )
    B = Symbolics.sparsejacobian(
        kkt_eqns[:],
        parameters,
    )

    syms_to_vals = merge(
        Dict(zip([xs; eq_duals; ineq_duals], [x_vals; eq_dual_vals[eq_dual_idxs]; ineq_dual_vals[ineq_dual_idxs]])),
        parameter_values,
    )

    # substitute in values
    Is, Js, Vs = findnz(A)
    vs = float.(Symbolics.value.(Symbolics.substitute(Vs, syms_to_vals)))
    a = sparse(Is, Js, vs, size(A)...)

    Is, Js, Vs = findnz(B)
    vs = float.(Symbolics.value.(Symbolics.substitute(Vs, syms_to_vals)))
    b = Array(sparse(Is, Js, vs, size(B)...)) # no sparse rhs solver, need to make dense

   rank(a)
    t = qr(a)

    -a \ b # sensitivities
end

export kkt
