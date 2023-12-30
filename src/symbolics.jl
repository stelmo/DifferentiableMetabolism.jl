
#=
Extend the substitute function from symbolics to convert parameters to numbers.
Converts all the ParameterXXX types to their XXX types normally found in
ConstraintTrees.
=#

Symbolics.substitute(x::ParameterLinearValue, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.LinearValue(
        x.idxs,
        Symbolics.value.(Symbolics.substitute(x.weights, rule)),
    )

Symbolics.substitute(x::ParameterBetween, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.Between(
        Symbolics.value.(Symbolics.substitute(x.lower, rule)),
        Symbolics.value.(Symbolics.substitute(x.upper, rule)),
    )

Symbolics.substitute(x::ParameterEqualTo, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.EqualTo(Symbolics.value.(Symbolics.substitute(x.equal_to, rule)))

Symbolics.substitute(x::ConstraintTrees.Between, rule::Dict{Symbolics.Num,Float64}) = x

Symbolics.substitute(x::ConstraintTrees.EqualTo, rule::Dict{Symbolics.Num,Float64}) = x

Symbolics.substitute(x::ConstraintTrees.LinearValue, rule::Dict{Symbolics.Num,Float64}) = x

Symbolics.substitute(x::Nothing, rule::Dict{Symbolics.Num,Float64}) = nothing

Symbolics.substitute(x::ConstraintTrees.Constraint, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.Constraint(
        Symbolics.substitute(ConstraintTrees.value(x), rule),
        Symbolics.substitute(ConstraintTrees.bound(x), rule),
    )

Symbolics.substitute(x::ConstraintTrees.ConstraintTree, rule::Dict{Symbolics.Num,Float64}) =
    ConstraintTrees.tree_map(
        x,
        c -> Symbolics.substitute(c, rule),
        ConstraintTrees.Constraint,
    )
