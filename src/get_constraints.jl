"""
$(TYPEDSIGNATURES)

Return all the equality constraints of `m` as a tuple `({Parameter}LinearValue,
value)` representing `{P}LV == value` for each entry.
"""
function equality_constraints(m::ConstraintTrees.ConstraintTree)
    sink = Vector{
        Tuple{Union{ParameterLinearValue,ConstraintTrees.LinearValue},Symbolics.Num},
    }()
    get_equality_constraints(m, sink)
    sink
end

export equality_constraints

function get_equality_constraints(m::ConstraintTrees.ConstraintTree, sink)
    get_equality_constraints.(values(m), Ref(sink))
end

function get_equality_constraints(c::ConstraintTrees.Constraint, sink)
    if c.bound isa ConstraintTrees.EqualTo || c.bound isa ParameterEqualTo
        push!(sink, (ConstraintTrees.value(c), Symbolics.Num(c.bound.equal_to)))
    end
end

"""
$(TYPEDSIGNATURES)

Return all the inequality constraints of `m` as a tuple `({Parameter}LinearValue,
lower, upper)` representing `lower ≤ {P}LV ≤ upper` for each entry.
"""
function inequality_constraints(m::ConstraintTrees.ConstraintTree)
    sink = Vector{
        Tuple{
            Union{ParameterLinearValue,ConstraintTrees.LinearValue},
            Symbolics.Num,
            Symbolics.Num,
        },
    }()
    get_inequality_constraints(m, sink)
    sink
end

export inequality_constraints

function get_inequality_constraints(m::ConstraintTrees.ConstraintTree, sink)
    get_inequality_constraints.(values(m), Ref(sink))
end

function get_inequality_constraints(c::ConstraintTrees.Constraint, sink)
    if c.bound isa ConstraintTrees.Between || c.bound isa ParameterBetween
        push!(
            sink,
            (
                ConstraintTrees.value(c),
                Symbolics.Num(c.bound.lower),
                Symbolics.Num(c.bound.upper),
            ),
        )
    end
end
