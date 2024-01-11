
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
=#


function remove_linearly_dependent_constraints(eqs, parameter_values, xs)

    #= 
    Filter out linearly dependent constraints using QR decomposition. Since the
    problem solved, assume there are no contradictory constraints. 

    See: https://math.stackexchange.com/questions/748500/how-to-find-linearly-independent-columns-in-a-matrix

    Make use of the fact that sparse QR returns a staircase profile with column
    ordering by default. The default tolerance of what is a 0 seems okay to rely
    on. Could be a source of bugs though...
    =#
    Is = Int[]
    Js = Int[]
    Vs = Float64[]
    for (i, eq) in enumerate(eqs)
        lhs = first(eq)
        append!(Is, fill(i, length(lhs.idxs)))
        append!(Js, lhs.idxs)
        append!(Vs, Symbolics.value.(Symbolics.substitute(lhs.weights, parameter_values)))
    end
    mat_sparse_transposed = sparse(Js, Is, Vs) # do transpose here for QR

    t = qr(mat_sparse_transposed)
    max_lin_indep_columns = 0
    for i in axes(t.R, 2) # depends on preordered QR!
        Is, _ = findnz(t.R[:, i])
        if maximum(Is) == i
            max_lin_indep_columns = i
        end
    end
    lin_indep_rows = t.pcol[1:max_lin_indep_columns] # undo permumation

    # final equality constraints
    [ConstraintTrees.substitute(lhs, xs) - rhs for (lhs, rhs) in eqs[lin_indep_rows]],
    lin_indep_rows
end

export remove_linearly_dependent_constraints

"""
$(TYPEDSIGNATURES)

Differentiate a model `m` with respect to `parameters` which take on values
`parameter_values` in the optimal solution with respect to the `objective`. The
primal variables `x_vals`, and the dual variable values `eq_dual_vals` and
`ineq_dual_vals` need to be supplied. 

Internally, primal variables with value `abs(x) ≤ primal_zero_tol` are removed
from the computation, and their sensitivities are not calculated.  
"""
function differentiate(
    m::ConstraintTrees.ConstraintTree,
    objective::ConstraintTrees.Value,
    x_vals::Vector{Float64},
    eq_dual_vals::Vector{Float64},
    ineq_dual_vals::Vector{Float64},
    parameter_values::Dict{Symbolics.Num,Float64},
    parameters::Vector{Symbolics.Num}; # might not diff wrt all params
    rank_zero_tol = 1e-8,
)
    # create symbolic values of the primal and dual variables
    Symbolics.@variables x[1:ConstraintTrees.var_count(m)]
    xs = collect(x) # to make overloads in DiffMet work correctly

    # objective
    f = ConstraintTrees.substitute(objective, xs)

    # equality constraints
    # E * x - b = H = 0
    H, lin_indep_rows =
        remove_linearly_dependent_constraints(equality_constraints(m), parameter_values, xs)

    # inequality constraints (must be built the same as in solver.jl)
    # M * x - h = G ≤ 0
    ineqs = inequality_constraints(m)
    G = [ConstraintTrees.substitute(lhs, xs) - rhs for (lhs, rhs) in ineqs]


    # creaty symbolic variables for the duals, but only those that are required
    Symbolics.@variables eq_duals[1:length(H)] ineq_duals[1:length(G)]

    #=
    Do all the manipulations manually. This is much faster than using the
    builtin functions.

    kkt_eqns = [ 
        ∇ₓfᵀ - ∇ₓHᵀ ν  - ∇ₓGᵀ λ # negatives because of KKT formulation in JuMP
        H
        G .* ineq_duals
    ]
    =#
    eq1 = Symbolics.sparsejacobian([f], xs)'

    Is, Js, Vs = findnz(Symbolics.sparsejacobian(H, xs))
    eq2 = zeros(Symbolics.Num, size(eq1, 1))
    for (i, j, v) in zip(Js, Is, Vs) # transpose
        eq2[i] += v * eq_duals[j]
    end

    Is, Js, Vs = findnz(Symbolics.sparsejacobian(G, xs))
    eq3 = zeros(Symbolics.Num, size(eq1, 1))
    for (i, j, v) in zip(Js, Is, Vs) # transpose
        eq3[i] += v * ineq_duals[j]
    end

    eq4 = zeros(Symbolics.Num, size(ineq_duals, 1))
    for i in eachindex(G) # transpose
        eq4[i] += G[i] * ineq_duals[i]
    end

    kkt_eqns = [ # negatives because of KKT formulation in JuMP
        eq1 - eq2 - eq3
        H
        eq4
    ]

    A = Symbolics.sparsejacobian(kkt_eqns[:], [xs; eq_duals; ineq_duals])
    B = Symbolics.sparsejacobian(kkt_eqns[:], parameters)

    # symbolic values at the optimal solution incl parameters
    syms_to_vals = merge(
        Dict(
            zip(
                [xs; eq_duals; ineq_duals],
                [x_vals; eq_dual_vals[lin_indep_rows]; ineq_dual_vals],
            ),
        ),
        parameter_values,
    )

    # substitute in values
    Is, Js, Vs = findnz(A)
    vs = float.(Symbolics.value.(Symbolics.substitute(Vs, syms_to_vals)))
    a = sparse(Is, Js, vs, size(A)...)
    # rank(a)

    # _,_,vs = findnz(a)
    # minimum(abs, vs)

    # size(eq1)
    # af = a[1:size(eq1, 1), :]
    # rank(af)

    # size(H)
    # ah = a[(1+size(eq1, 1)):(size(H, 1)+size(eq1, 1)), :]
    # rank(ah)

    # size(eq4)
    # ag = a[(1+size(eq1, 1)+size(H, 1)):end, :]
    # rank(ag)

    # t = qr(a)
    # max_lin_indep_columns = 0
    # for i in axes(t.R, 2) # depends on preordered QR!
    #     Is, _ = findnz(t.R[:, i])
    #     if maximum(Is) == i
    #         max_lin_indep_columns = i
    #     end
    # end
    # lin_indep_rows = t.pcol[1:max_lin_indep_columns] # undo permumation
    # setdiff(collect(1:size(a, 2)), lin_indep_rows)


    Is, Js, Vs = findnz(B)
    vs = float.(Symbolics.value.(Symbolics.substitute(Vs, syms_to_vals)))
    b = Array(sparse(Is, Js, vs, size(B)...)) # no sparse rhs solver, need to make dense

    # calculate sensitivities
    if size(a, 1) >= size(a, 2)
        #=
        If a is rectangular (more equations than variables), then this should
        still work because the equations should not be in conflict (in an ideal
        world).
        =#
        if rank(a; tol = rank_zero_tol) < size(a, 2)
            throw(
                ArgumentError(
                    "The optimal solution is not unique, the derivatives cannot be calculated.",
                ),
            )
        else
            -a \ b
        end

    else # something went wrong
        zeros(0, 0)
    end
end

export differentiate
