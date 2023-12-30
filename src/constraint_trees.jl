
#=
Extend various functions in ConstraintTrees to work with their Parameter based equivalents.
=#

ConstraintTrees.incr_var_idxs(x::ParameterLinearValue, incr::Int) = ParameterLinearValue(
    idxs = ConstraintTrees.incr_var_idx.(x.idxs, incr),
    weights = x.weights,
)

ConstraintTrees.var_count(x::ParameterLinearValue) = isempty(x.idxs) ? 0 : last(x.idxs)

# substitute in the variables as symbolic numbers - useful to construct the KKT function
ConstraintTrees.substitute(x::ParameterLinearValue, y::Vector{Symbolics.Num}) = sum(
    (idx == 0 ? x.weights[i] : x.weights[i] * y[idx] for (i, idx) in enumerate(x.idxs)),
    init = Symbolics.Num(0.0),
)

ConstraintTrees.substitute(x::ParameterQuadraticValue, y::Vector{Symbolics.Num}) = sum(
    (
        let (idx1, idx2) = x.idxs[i]
            (idx1 == 0 ? 1.0 : y[idx1]) * (idx2 == 0 ? 1.0 : y[idx2]) * w
        end for (i, w) in enumerate(x.weights)
    ),
    init = Symbolics.Num(0.0),
)
