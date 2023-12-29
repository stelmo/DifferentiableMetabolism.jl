C.incr_var_idxs(x::ParameterLinearValue, incr::Int) = ParameterLinearValue(idxs = C.incr_var_idx.(x.idxs, incr), weights = x.weights)

C.var_count(x::ParameterLinearValue) = isempty(x.idxs) ? 0 : last(x.idxs)

# substitute in the variables (to differentiate in KKT)
C.substitute(x::ParameterLinearValue, y::Vector{S.Num}) = sum(
    (idx == 0 ? x.weights[i] : x.weights[i] * y[idx] for (i, idx) in enumerate(x.idxs)),
    init = S.Num(0.0),
)

C.constraint_values(x::C.Constraint, y::Vector{S.Num}) = C.substitute(C.value(x), y)

C.constraint_values(x::C.ConstraintTree, y::Vector{S.Num}) = C.tree_map(x, c -> C.substitute(C.value(c), y), S.Num)

# substitute in ALL the parameters (to solve normally)
C.substitute(x::ParameterLinearValue, rule::Dict{S.Num, Float64}) = C.LinearValue(x.idxs, S.value.(S.substitute(x.weights, rule)))

C.substitute(x::ParameterBetween, rule::Dict{S.Num, Float64}) = C.Between(S.value.(S.substitute(x.lower, rule)), S.value.(S.substitute(x.upper, rule)))

# these don't have parameters, nothing needs to happen
C.substitute(x::C.Between, rule::Dict{S.Num, Float64}) = x
C.substitute(x::C.EqualTo, rule::Dict{S.Num, Float64}) = x
C.substitute(x::C.LinearValue, rule::Dict{S.Num, Float64}) = x
C.substitute(x::Nothing, rule::Dict{S.Num, Float64}) = nothing

C.substitute(x::C.Constraint, rule::Dict{S.Num, Float64}) = C.Constraint(
    C.substitute(C.value(x), rule),
    C.substitute(C.bound(x), rule),
)

C.substitute(x::C.ConstraintTree, rule::Dict{S.Num, Float64}) = C.tree_map(x, c -> C.substitute(c, rule), Any) # TODO what should the type be?
