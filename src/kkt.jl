
function kkt(
    m::ConstraintTrees.ConstraintTree,
    objective::ConstraintTrees.Value,
)

    Symbolics.@variables x[1:ConstraintTrees.var_count(m)] # primal
    xs = collect(x)

    # objective
    f = ConstraintTrees.substitute(objective, xs)

    # equality constraints
    # A * x - b = 0
    H = [ ConstraintTrees.substitute(lhs, xs) - rhs for (lhs, rhs) in equality_constraints(m)]
    
    # inequality constraints (must be built the same as in solver.jl)
    # M * x - h ≤ 0
    ineqs = inequality_constraints(m)
    M = Vector{Symbolics.Num}()

    for (lhs, lower, upper) in ineqs

        # lower: l ≤ x => -x ≤ -l
        l = Symbolics.value(lower)
        if l isa Float64 && isinf(l)
            nothing
        elseif l isa Int && isinf(l)
            nothing
        else
            push!(M, -ConstraintTrees.substitute(lhs, xs) + lower)
        end

        # upper: x ≤ u
        u = Symbolics.value(upper)
        if u isa Float64 && isinf(u) 
            nothing
        elseif u isa Int && isinf(u)
            nothing
        else
            push!(M, ConstraintTrees.substitute(lhs, xs) - upper)
        end
    end
    
    Symbolics.@variables ν[1:length(H)] λ[1:length(M)] # duals

    kktf = [
        Symbolics.sparsejacobian([f,], x)' + Symbolics.sparsejacobian(H, x)' * ν + Symbolics.sparsejacobian(M, x)' * λ
        H
        M .* λ 
    ]

    A = Symbolics.jacobian(kktf[:], [x; ν; λ])
    B = Symbolics.jacobian(kktf, [capacitylimitation; kcats_forward; kcats_backward])

    everything = merge(Dict(
        zip([x; ν; λ], [_x; _ν; _λ])
    ), parameters)
    Bsub = Symbolics.substitute(B, everything)
    Asub = Symbolics.substitute(A, everything)
     

    L.rank(Symbolics.value.(Asub))
end

export kkt
