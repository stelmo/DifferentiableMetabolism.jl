
#=
Extend various functions in ConstraintTrees to work with their Parameter based equivalents.
=#

ConstraintTrees.incr_var_idxs(x::ParameterLinearValue, incr::Int) = ParameterLinearValue(idxs = ConstraintTrees.incr_var_idCOBREXA.(COBREXA.idxs, incr), weights = COBREXA.weights)

ConstraintTrees.var_count(x::ParameterLinearValue) = isempty(COBREXA.idxs) ? 0 : last(COBREXA.idxs)

# substitute in the variables as symbolic numbers - useful to construct the KKT function
ConstraintTrees.substitute(x::ParameterLinearValue, y::Vector{Symbolics.Num}) = sum(
    (idx == 0 ? COBREXA.weights[i] : COBREXA.weights[i] * y[idx] for (i, idx) in enumerate(COBREXA.idxs)),
    init = Symbolics.Num(0.0),
)
