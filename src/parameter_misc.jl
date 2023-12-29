C.incr_var_idxs(x::ParameterLinearValue, incr::Int) =
    ParameterLinearValue(idxs = C.incr_var_idx.(x.idxs, incr), weights = x.weights)

C.var_count(x::ParameterLinearValue) = isempty(x.idxs) ? 0 : last(x.idxs)

# substitute in ALL the parameters (to solve normally)
S.substitute(x::ParameterLinearValue, rule::Dict{S.Num,Float64}) =
    C.LinearValue(x.idxs, S.value.(S.substitute(x.weights, rule)))

S.substitute(x::ParameterBetween, rule::Dict{S.Num,Float64}) =
    C.Between(S.value.(S.substitute(x.lower, rule)), S.value.(S.substitute(x.upper, rule)))

S.substitute(x::ParameterEqualTo, rule::Dict{S.Num,Float64}) =
    C.EqualTo(S.value.(S.substitute(x.equal_to, rule)))

# these don't have parameters, nothing needs to happen
S.substitute(x::C.Between, rule::Dict{S.Num,Float64}) = x
S.substitute(x::C.EqualTo, rule::Dict{S.Num,Float64}) = x
S.substitute(x::C.LinearValue, rule::Dict{S.Num,Float64}) = x
S.substitute(x::Nothing, rule::Dict{S.Num,Float64}) = nothing

S.substitute(x::C.Constraint, rule::Dict{S.Num,Float64}) =
    C.Constraint(S.substitute(C.value(x), rule), S.substitute(C.bound(x), rule))

S.substitute(x::C.ConstraintTree, rule::Dict{S.Num,Float64}) =
    C.tree_map(x, c -> S.substitute(c, rule), C.Constraint)


# substitute in the variables (to differentiate in KKT)
C.substitute(x::ParameterLinearValue, y::Vector{S.Num}) = sum(
    (idx == 0 ? x.weights[i] : x.weights[i] * y[idx] for (i, idx) in enumerate(x.idxs)),
    init = S.Num(0.0),
)

function equality_constraints(c::C.ConstraintTree)
    sink = Vector{Tuple{Union{ParameterLinearValue,C.LinearValue},S.Num}}()
    get_equality_constraints(c, sink)
    sink
end

export equality_constraints

function get_equality_constraints(c::C.ConstraintTree, sink)
    get_equality_constraints.(values(c), Ref(sink))
end

function get_equality_constraints(c::C.Constraint, sink)
    if c.bound isa C.EqualTo || c.bound isa ParameterEqualTo
        push!(sink, (C.value(c), S.Num(c.bound.equal_to)))
    end
end

function inequality_constraints(c::C.ConstraintTree)
    sink = Vector{Tuple{Union{ParameterLinearValue,C.LinearValue},S.Num,S.Num}}()
    get_inequality_constraints(c, sink)
    sink
end

export inequality_constraints

function get_inequality_constraints(c::C.ConstraintTree, sink)
    get_inequality_constraints.(values(c), Ref(sink))
end

function get_inequality_constraints(c::C.Constraint, sink)
    if c.bound isa C.Between || c.bound isa ParameterBetween
        push!(sink, (C.value(c), S.Num(c.bound.lower), S.Num(c.bound.upper)))
    end
end
