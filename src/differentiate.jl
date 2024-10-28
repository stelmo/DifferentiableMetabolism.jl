
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

"""
$(TYPEDSIGNATURES)

TODO
"""
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

export differentiate_prepare_kkt

"""
$(TYPEDSIGNATURES)

Prepare a model `m` with `objective` for differentiation with respect to `parameters`.
"""
function differentiate_prepare_kkt(
    m::ConstraintTrees.ConstraintTree,
    objective::ConstraintTrees.Value,
    parameters::Vector{Symbol}; # might not diff wrt all params
)
    # create symbolic values of the primal and dual variables
    xs = FastDifferentiation.make_variables(:x, ConstraintTrees.var_count(m))

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
    eq_duals = FastDifferentiation.make_variables(:eq_duals, length(H))
    ineq_duals = FastDifferentiation.make_variables(:ineq_duals, length(G))

    #=
    Do all the manipulations manually. This is much faster than using the
    builtin functions.

    kkt_eqns = [ 
        ∇ₓfᵀ - ∇ₓHᵀ ν  - ∇ₓGᵀ λ # negatives because of KKT formulation in JuMP
        H
        G .* ineq_duals
    ]

    Note, make sure all the lazy operations are expanded to avoid running into bugs:
    https://github.com/JuliaFastDifferentiation/FastDifferentiation.jl/issues/518
    https://github.com/JuliaFastDifferentiation/FastDifferentiation.jl/issues/498
    =#
    eq1 = FastDifferentiation.jacobian([f], xs)[1,:]

    Is, Js, Vs = SparseArrays.findnz(FastDifferentiation.sparse_jacobian(H, xs))
    eq2 = zeros(Expression, size(eq1, 1))
    for (i, j, v) in zip(Js, Is, Vs) # transpose
        eq2[i] += v * eq_duals[j]
    end

    Is, Js, Vs = SparseArrays.findnz(FastDifferentiation.sparse_jacobian(G, xs))
    eq3 = zeros(Expression, size(eq1, 1))
    for (i, j, v) in zip(Js, Is, Vs) # transpose
        eq3[i] += v * ineq_duals[j]
    end

    eq4 = zeros(Expression, size(ineq_duals, 1))
    for i in eachindex(G) # transpose
        eq4[i] += G[i] * ineq_duals[i]
    end

    kkt_eqns = [ # negatives because of KKT formulation in JuMP
        eq1 - eq2 - eq3
        H
        eq4
    ]

    A = FastDifferentiation.sparse_jacobian(kkt_eqns, [xs; eq_duals; ineq_duals])
    B = FastDifferentiation.sparse_jacobian(kkt_eqns, FastDifferentiation.Node.(parameters))

    return (A, B, xs, eq_duals, ineq_duals, parameters), variable_order(m)
end

export differentiate_prepare_kkt

function differentiate_solution(
    (A, B, xs, eq_duals, ineq_duals, parameters),
    x_vals::Vector{Float64},
    eq_dual_vals::Vector{Float64},
    ineq_dual_vals::Vector{Float64},
    parameter_values::Dict{Symbol,Float64};
    scale = false, # scale sensitivities
)

    # symbolic values at the optimal solution incl parameters
    syms_to_vals = merge(
        Dict(zip((x.node_value for x in [xs; eq_duals; ineq_duals]), [x_vals; eq_dual_vals; ineq_dual_vals])),
        parameter_values,
    )

    # substitute in values
    Is, Js, Vs = SparseArrays.findnz(A)
    vs = float.(substitute.(Vs, Ref(k -> syms_to_vals[k])))
    a = SparseArrays.sparse(Is, Js, vs, size(A)...)
    indep_rows = findall_indeps_qr(a) # find independent rows, prevent singularity issues with \
    a_indep = a[indep_rows, :]

    #=
    If a is rectangular (more equations than variables), then the above should
    be sufficient, because the equations should not be in conflict (in an ideal
    world). 
    =#

    Is, Js, Vs = SparseArrays.findnz(B)
    vs = float.(substitute.(Vs, Ref(k -> syms_to_vals[k])))
    b = Array(SparseArrays.sparse(Is, Js, vs, size(B)...)) # no sparse rhs solver, need to make dense
    b_indep = b[indep_rows, :]

    c = -a_indep \ b_indep # sensitivities, unscaled

    # get primal variable sensitivities only
    if scale
        ([parameter_values[p] for p in parameters]' .* c[1:length(xs), :] ./ x_vals)
    else
        c[1:length(xs),:]
    end
end

export differentiate_solution

"""
$(TYPEDSIGNATURES)

TODO
"""
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
