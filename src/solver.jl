
function optimization_model_with_parameters(
    m::C.ConstraintTree,
    parameters::Dict{S.Num,Float64};
    objective::C.Value,
    optimizer,
    sense,
)
    model = J.Model(optimizer)

    J.@variable(model, x[1:C.var_count(m)])
    J.@objective(model, sense, C.substitute(S.substitute(objective, parameters), x))

    eqs = equality_constraints(m)
    ineqs = inequality_constraints(m)

    # A * x = b
    Ib = Int64[]
    Vb = Float64[]
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    for (i, (val, rhs)) in enumerate(eqs)
        rhs = S.substitute(rhs, parameters)
        a = S.substitute(val, parameters)

        if S.value(rhs) != 0
            push!(Ib, i)
            push!(Vb, S.value(rhs))
        end
        append!(Is, fill(i, length(a.idxs)))
        append!(Js, a.idxs)
        append!(Vs, a.weights)
    end
    A = sparse(Is, Js, Vs, length(eqs), C.var_count(m))
    b = sparsevec(Ib, Vb, length(eqs))
    J.@constraint(model, eqcons, A * x .== b)

    # M * x ≤ h
    Ih = Int64[]
    Vh = Float64[]
    Is = Int64[]
    Js = Int64[]
    Vs = Float64[]
    k = 0
    for (val, lower, upper) in ineqs
        lower = S.substitute(lower, parameters)
        upper = S.substitute(upper, parameters)
        a = S.substitute(val, parameters)

        # lower: l ≤ x => -x ≤ -l
        if !isinf(S.value(lower))
            k += 1
            if S.value(lower) != 0
                push!(Ih, k)
                push!(Vh, -S.value(lower))
            end
            append!(Is, fill(k, length(a.idxs)))
            append!(Js, a.idxs)
            append!(Vs, -a.weights)
        end

        # upper: x ≤ u
        if !isinf(S.value(upper))
            k += 1
            if S.value(upper) != 0
                push!(Ih, k)
                push!(Vh, S.value(upper))
            end
            append!(Is, fill(k, length(a.idxs)))
            append!(Js, a.idxs)
            append!(Vs, a.weights)
        end
    end
    M = sparse(Is, Js, Vs, k, C.var_count(m))
    h = sparsevec(Ih, Vh, k)
    J.@constraint(model, ineqcons, M * x .<= h)

    return model
end

export optimization_model_with_parameters

function optimized_constraints_with_parameters(
    m::C.ConstraintTree,
    parameters::Dict{S.Num,Float64};
    modifications = [],
    objective::C.Value,
    optimizer,
    sense = X.Maximal,
    duals = false,
)
    om = optimization_model_with_parameters(m, parameters; objective, optimizer, sense)
    for m in modifications
        m(om)
    end
    J.optimize!(om)

    X.is_solved(om) ?
    (
        duals ?
        (
            C.constraint_values(S.substitute(m, parameters), J.value.(om[:x])),
            J.dual.(om[:eqcons]),
            J.dual.(om[:ineqcons]),
        ) : C.constraint_values(S.substitute(m, parameters), J.value.(om[:x]))
    ) : nothing
end

export optimized_constraints_with_parameters
