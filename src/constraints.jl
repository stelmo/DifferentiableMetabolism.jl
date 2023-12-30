
function equality_constraints(c::ConstraintTrees.ConstraintTree)
    sink = Vector{
        Tuple{Union{ParameterLinearValue,ConstraintTrees.LinearValue},Symbolics.Num},
    }()
    get_equality_constraints(c, sink)
    sink
end

export equality_constraints

function get_equality_constraints(c::ConstraintTrees.ConstraintTree, sink)
    get_equality_constraints.(values(c), Ref(sink))
end

function get_equality_constraints(c::ConstraintTrees.Constraint, sink)
    if c.bound isa ConstraintTrees.EqualTo || c.bound isa ParameterEqualTo
        push!(sink, (ConstraintTrees.value(c), Symbolics.Num(c.bound.equal_to)))
    end
end

function inequality_constraints(c::ConstraintTrees.ConstraintTree)
    sink = Vector{
        Tuple{
            Union{ParameterLinearValue,ConstraintTrees.LinearValue},
            Symbolics.Num,
            Symbolics.Num,
        },
    }()
    get_inequality_constraints(c, sink)
    sink
end

export inequality_constraints

function get_inequality_constraints(c::ConstraintTrees.ConstraintTree, sink)
    get_inequality_constraints.(values(c), Ref(sink))
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
