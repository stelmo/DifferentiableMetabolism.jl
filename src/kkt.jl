
function kkt(
    m::C.ConstraintTree,
    objective::C.Value,
)

    S.@variables x[1:C.var_count(m)] # primal
    xs = collect(x)

    # objective
    f = C.substitute(objective, xs)

    # equality constraints
    # A * x - b = 0
    H = [ C.substitute(lhs, xs) - rhs for (lhs, rhs) in equality_constraints(m)]
    
    # inequality constraints (must be built the same as in solver.jl)
    # M * x - h ≤ 0
    ineqs = inequality_constraints(m)
    M = Vector{S.Num}()

    for (lhs, lower, upper) in ineqs

        # lower: l ≤ x => -x ≤ -l
        l = S.value(lower)
        if l isa Float64 && isinf(l)
            nothing
        elseif l isa Int && isinf(l)
            nothing
        else
            push!(M, -C.substitute(lhs, xs) + lower)
        end

        # upper: x ≤ u
        u = S.value(upper)
        if u isa Float64 && isinf(u) 
            nothing
        elseif u isa Int && isinf(u)
            nothing
        else
            push!(M, C.substitute(lhs, xs) - upper)
        end
    end
    
    S.@variables ν[1:length(H)] λ[1:length(M)] # duals

    kktf = [
        S.sparsejacobian([f,], x)' + S.sparsejacobian(H, x)' * ν + S.sparsejacobian(M, x)' * λ
        H
        M .* λ 
    ]

    A = S.jacobian(kktf[:], [x; ν; λ])
    B = S.jacobian(kktf, [capacitylimitation; kcats_forward; kcats_backward])

    everything = merge(Dict(
        zip([x; ν; λ], [_x; _ν; _λ])
    ), parameters)
    Bsub = S.substitute(B, everything)
    Asub = S.substitute(A, everything)
     

    L.rank(S.value.(Asub))
end

export kkt
