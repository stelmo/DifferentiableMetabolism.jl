
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

Return all linearly dependent constraints in  `A`, using the QR decomposition.
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
    Is, Js, Vs = SA.findnz(A)

    a = SA.sparse(Js, Is, Vs)


    t = LA.qr(a)  # do transpose here for QR
    max_lin_indep_columns = 0
    for i in axes(t.R, 2) # depends on preordered QR!
        Is, _ = SA.findnz(t.R[:, i])
        if isempty(Is) || maximum(Is) != i
            break
        else
            max_lin_indep_columns = i
        end
    end

    t.pcol[1:max_lin_indep_columns] # undo permutation
end

"""
$(TYPEDSIGNATURES)

Prepare a model `m` with `objective` for differentiation with respect to
`parameters`. 

This is the most time consuming aspect of differentiation. It pays off to do
this separately  if the same model will be differentiated multiple times.
"""
function differentiate_prepare_kkt(
    m::C.ConstraintTree,
    objective::C.Value,
    parameters::Vector{Symbol}; # might not diff wrt all params
)
    # create symbolic values of the primal and dual variables
    primals = F.make_variables(:x, C.var_count(m))

    # objective
    f = C.substitute(objective, primals)

    # equality constraints
    # E * x - b = H = 0
    eqs = equality_constraints(m)
    H = [C.substitute(lhs, primals) - rhs for (lhs, rhs) in eqs]

    # inequality constraints (must be built the same as in solver.jl)
    # M * x - h = G ≤ 0
    ineqs = inequality_constraints(m)
    G = [C.substitute(lhs, primals) - rhs for (lhs, rhs) in ineqs]

    # creaty symbolic variables for the duals, but only those that are required
    eq_duals = F.make_variables(:eq_duals, length(H))
    ineq_duals = F.make_variables(:ineq_duals, length(G))

    #=
    Do all the manipulations manually. This is much faster than using the
    builtin functions.

    kkt_eqns = [
        ∇ₓfᵀ - ∇ₓHᵀ ν  - ∇ₓGᵀ λ # negatives because of KKT formulation in JuMP
        H
        G .* ineq_duals
    ]

    Note, make sure all the lazy operations are expanded to avoid running into bugs...
    =#
    eq1 = F.jacobian([f], primals)[1, :]

    Is, Js, Vs = SA.findnz(F.sparse_jacobian(H, primals))
    eq2 = zeros(Ex, size(eq1, 1))
    for (i, j, v) in zip(Js, Is, Vs) # transpose
        eq2[i] += v * eq_duals[j]
    end

    Is, Js, Vs = SA.findnz(F.sparse_jacobian(G, primals))
    eq3 = zeros(Ex, size(eq1, 1))
    for (i, j, v) in zip(Js, Is, Vs) # transpose
        eq3[i] += v * ineq_duals[j]
    end

    eq4 = zeros(Ex, size(ineq_duals, 1))
    for i in eachindex(G) # transpose
        eq4[i] += G[i] * ineq_duals[i]
    end

    kkt_eqns = [ # negatives because of KKT formulation in JuMP
        eq1 - eq2 - eq3
        H
        eq4
    ]

    A = F.sparse_jacobian(kkt_eqns, [primals; eq_duals; ineq_duals])
    B = F.sparse_jacobian(kkt_eqns, F.Node.(parameters))

    return (A, B, primals, eq_duals, ineq_duals, parameters), variable_order(m)
end

"""
$(TYPEDSIGNATURES)

Differentiate a solution of a model. The first argument is the output of [`differentiate_prepare_kkt`](@ref), and is a tuple of the deconstructed model.
The following arguments (`primal_vals`, `eq_dual_vals`, `ineq_dual_vals`) are outputs of [`optimized_constraints_with_parameters`](@ref).
`parameter_values`
"""
function differentiate_solution(
    (A, B, primals, eq_duals, ineq_duals, parameters),
    primal_vals::Vector{Float64},
    eq_dual_vals::Vector{Float64},
    ineq_dual_vals::Vector{Float64},
    parameter_values::Dict{Symbol,Float64};
    scale = false, # scale sensitivities
)

    # symbolic values at the optimal solution incl parameters
    syms_to_vals = merge(
        Dict(
            zip(
                (x.node_value for x in [primals; eq_duals; ineq_duals]),
                [primal_vals; eq_dual_vals; ineq_dual_vals],
            ),
        ),
        parameter_values,
    )

    # substitute in values
    Is, Js, Vs = SA.findnz(A)
    vs = float.(substitute.(Vs, Ref(k -> syms_to_vals[k])))
    a = SA.sparse(Is, Js, vs, size(A)...)
    indep_rows = findall_indeps_qr(a) # find independent rows, prevent singularity issues with \
    a_indep = a[indep_rows, :]

    #=
    If a is rectangular (more equations than variables), then the above should
    be sufficient, because the equations should not be in conflict (in an ideal
    world).
    =#

    Is, Js, Vs = SA.findnz(B)
    vs = float.(substitute.(Vs, Ref(k -> syms_to_vals[k])))
    b = Array(SA.sparse(Is, Js, vs, size(B)...)) # no sparse rhs solver, need to make dense
    b_indep = b[indep_rows, :]

    c = -a_indep \ b_indep # sensitivities, unscaled

    # get primal variable sensitivities only
    if scale
        (
            [parameter_values[p] for p in parameters]' .* c[1:length(primals), :] ./ primal_vals
        )
    else
        c[1:length(primals), :]
    end
end

"""
$(TYPEDSIGNATURES)

Return the names of variables in `m`.

NOTE: this function assumes that all variables are bounded explicitly in the
model. `nothing` bounds are ignored.
"""
function variable_order(m)
    c = []
    ff(p, x::C.ConstraintTree) = nothing
    ff(p, x::C.Constraint) = begin
        if length(x.value.idxs) == 1 && !isnothing(x.bound) #TODO assumes that all variables are bounded somehow!!
            push!(c, (first(x.value.idxs), p))
        end
    end

    C.itraverse(ff, m)

    idxs = first.(c)
    _idxs = sortperm(idxs)
    last.(c)[_idxs]
end
