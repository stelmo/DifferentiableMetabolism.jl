
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

function findall_indeps_qr(A)
    #= 
    Filter out linearly dependent constraints using QR decomposition. Since the
    problem solved, assume there are no contradictory constraints. 

    See: https://math.stackexchange.com/questions/748500/how-to-find-linearly-independent-columns-in-a-matrix

    Make use of the fact that sparse QR returns a staircase profile with column
    ordering by default. The default tolerance of what is a 0 seems okay to rely
    on. Could be a source of bugs though...
    =#
    Is, Js, Vs = SparseArrays.findnz(A)

    a = SparseArrays.sparse(Js, Is, Vs)


    t = LinearAlgebra.qr(a)  # do transpose here for QR
    max_lin_indep_columns = 0
    for i in axes(t.R, 2) # depends on preordered QR!
        Is, _ = SparseArrays.findnz(t.R[:, i])
        if isempty(Is) || maximum(Is) != i
            break
        else
            max_lin_indep_columns = i
        end
    end

    t.pcol[1:max_lin_indep_columns] # undo permumation
end

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
    scale = false, # scale sensitivities
)
    # create symbolic values of the primal and dual variables
    Symbolics.@variables x[1:ConstraintTrees.var_count(m)]
    xs = collect(x) # to make overloads in DiffMet work correctly

    # objective
    f = ConstraintTrees.substitute(objective, xs)

    # equality constraints
    # E * x - b = H = 0
    eqs = equality_constraints(m)
    H = [ConstraintTrees.substitute(lhs, xs) - rhs for (lhs, rhs) in eqs]
    # H, lin_indep_rows = remove_linearly_dependent_constraints2(eqs, parameter_values, xs)

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

    Note, make sure all the lazy operations are expanded to avoid running into bugs:
    https://github.com/JuliaSymbolics/Symbolics.jl/issues/518
    https://github.com/JuliaSymbolics/Symbolics.jl/issues/498
    =#
    eq1 = Symbolics.sparsejacobian([f], xs)'

    Is, Js, Vs = SparseArrays.findnz(Symbolics.sparsejacobian(H, xs))
    eq2 = zeros(Symbolics.Num, size(eq1, 1))
    for (i, j, v) in zip(Js, Is, Vs) # transpose
        eq2[i] += v * eq_duals[j]
    end

    Is, Js, Vs = SparseArrays.findnz(Symbolics.sparsejacobian(G, xs))
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
        Dict(zip([xs; eq_duals; ineq_duals], [x_vals; eq_dual_vals; ineq_dual_vals])),
        parameter_values,
    )

    # substitute in values
    Is, Js, Vs = SparseArrays.findnz(A)
    vs = float.(Symbolics.value.(Symbolics.substitute(Vs, syms_to_vals)))
    a = SparseArrays.sparse(Is, Js, vs, size(A)...)
    indep_rows = findall_indeps_qr(a) # find independent rows, prevent singularity issues with \
    a_indep = a[indep_rows, :]

    #=
    If a is rectangular (more equations than variables), then the above should
    be sufficient, because the equations should not be in conflict (in an ideal
    world). 
    =#

    Is, Js, Vs = SparseArrays.findnz(B)
    vs = float.(Symbolics.value.(Symbolics.substitute(Vs, syms_to_vals)))
    b = Array(SparseArrays.sparse(Is, Js, vs, size(B)...)) # no sparse rhs solver, need to make dense
    b_indep = b[indep_rows, :]

    c = -a_indep \ b_indep # sensitivities, unscaled

    # get primal variable sensitivities only
    if scale
        sc = similar(c[1:length(xs), :])
        for i in axes(sc, 1)
            for j in axes(sc, 2)
                sc[i, j] =
                    c[1:length(xs), :][i, j] * (parameter_values[parameters[j]] / x_vals[i])
            end
        end

        sc, variable_order(m)
    else
        c[1:length(xs), :], variable_order(m)
    end
end

export differentiate

function variable_order(m)
    c = []
    ff(p, x::ConstraintTrees.ConstraintTree) = nothing
    ff(p, x::ConstraintTrees.Constraint) = begin
        if length(x.value.idxs) == 1 && !isnothing(x.bound) #TODO assumes that all variables are bounded somehow!!
            push!(c, (first(x.value.idxs), p))
        end
    end

    ConstraintTrees.itraverse(ff, m)

    idxs = first.(c)
    _idxs = sortperm(idxs)
    last.(c)[_idxs]
end

