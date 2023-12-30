
function kkt(m::ConstraintTrees.ConstraintTree, objective::ConstraintTrees.Value)

    Symbolics.@variables x[1:ConstraintTrees.var_count(m)] # primal
    xs = collect(x)

    # objective
    f = ConstraintTrees.substitute(objective, xs)

    # equality constraints
    # E * x - b = H = 0 
    H = [
        ConstraintTrees.substitute(lhs, xs) - rhs for (lhs, rhs) in equality_constraints(m)
    ]

    # inequality constraints (must be built the same as in solver.jl)
    # M * x - h = G ≤ 0
    ineqs = inequality_constraints(m)
    G = Vector{Symbolics.Num}()

    for (lhs, lower, upper) in ineqs

        # lower: l ≤ x => -x ≤ -l
        l = Symbolics.value(lower)
        if l isa Float64 && isinf(l)
            nothing
        else
            push!(G, -ConstraintTrees.substitute(lhs, xs) + lower)
        end

        # upper: x ≤ u
        u = Symbolics.value(upper)
        if u isa Float64 && isinf(u)
            nothing
        else
            push!(G, ConstraintTrees.substitute(lhs, xs) - upper)
        end
    end

    Symbolics.@variables eq_duals[1:length(H)] ineq_duals[1:length(G)]

    kkt_eqns = [
        Symbolics.sparsejacobian([f], x)' +
        Symbolics.sparsejacobian(H, x)' * eq_duals +
        Symbolics.sparsejacobian(G, x)' * ineq_duals
        H
        G .* ineq_duals
    ]

    A = Symbolics.sparsejacobian(kkt_eqns[:], [x; eq_duals; ineq_duals])
    B = Symbolics.sparsejacobian(kkt_eqns[:], [capacitylimitation; kcats_forward; kcats_backward])

    everything = merge(Dict(zip([x; eq_duals; ineq_duals], [_x; _eq_duals; _ineq_duals])), parameters)

    # substitute in
    Is, Js, Vs = findnz(A)
    _Asub = sparse(Is, Js, Symbolics.substitute(Vs, everything), size(A)...)
    Is, Js, Vs = findnz(Asub)
    vs = float.(Symbolics.value.(Vs))
    a = sparse(Is, Js, vs, size(_Asub)...)
    
    Is, Js, Vs = findnz(B)
    _Bsub = sparse(Is, Js, Symbolics.substitute(Vs, everything), size(B)...)
    Is, Js, Vs = findnz(_Bsub)
    vs = float.(Symbolics.value.(Vs))
    b = Array(sparse(Is, Js, vs, size(_Bsub)...))
    
    a\b
end

export kkt
