
#=
Copyright (c) 2024, Heinrich-Heine University Duesseldorf
Copyright (c) 2024, University of Luxembourg

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

"""
$(TYPEDSIGNATURES)

TODO
"""
constrained(x) = begin
    y = FastDifferentiation.value(x)
    y isa Float64 && isinf(y) && return false
    return true
end

"""
$(TYPEDSIGNATURES)

Return all the equality constraints of `m` as a tuple `({Parameter}LinearValue,
value)` representing `{P}LV == value` for each entry.
"""
function equality_constraints(m::ConstraintTrees.ConstraintTree)
    sink =
        Vector{Tuple{Union{ParameterLinearValue,ConstraintTrees.LinearValue},Expression}}()
    get_equality_constraints(m, sink)
    sink
end

"""
$(TYPEDSIGNATURES)

TODO
"""
function get_equality_constraints(m::ConstraintTrees.ConstraintTree, sink)
    get_equality_constraints.(values(m), Ref(sink))
end

"""
$(TYPEDSIGNATURES)

TODO
"""
function get_equality_constraints(c::ConstraintTrees.Constraint, sink)
    if c.bound isa ConstraintTrees.EqualTo || c.bound isa ParameterEqualTo
        push!(sink, (ConstraintTrees.value(c), Expression(c.bound.equal_to)))
    end
end

"""
$(TYPEDSIGNATURES)

Return all the inequality constraints of `m` as a tuple of bounds converted to
the form `({Parameter}LinearValue, upper)` representing `{P}LV ≤ upper` for each
entry.
"""
function inequality_constraints(m::ConstraintTrees.ConstraintTree)
    sink = Vector{
        Tuple{
            Union{ParameterLinearValue,ConstraintTrees.LinearValue},
            Expression,
            Expression,
        },
    }()
    get_inequality_constraints(m, sink)

    [
        [ # lower bounds to l ≤ x → -x ≤ -l
            (-val, -lower) for (val, lower, _) in sink if constrained(lower)
        ]
        [ # upper bounds are already in the correct format
            (val, upper) for (val, _, upper) in sink if constrained(upper)
        ]
    ]
end

"""
$(TYPEDSIGNATURES)

TODO
"""
function get_inequality_constraints(m::ConstraintTrees.ConstraintTree, sink)
    get_inequality_constraints.(values(m), Ref(sink))
end

"""
$(TYPEDSIGNATURES)

TODO
"""
function get_inequality_constraints(c::ConstraintTrees.Constraint, sink)
    if c.bound isa ConstraintTrees.Between || c.bound isa ParameterBetween
        push!(
            sink,
            (
                ConstraintTrees.value(c),
                Expression(c.bound.lower),
                Expression(c.bound.upper),
            ),
        )
    end
end
