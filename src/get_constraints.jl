
#=
Copyright (c) 2025, Heinrich-Heine University Duesseldorf
Copyright (c) 2025, University of Luxembourg

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

Check if the argument is constrained (i.e. bound not infinity) or not.
"""
is_constrained(x) = begin
    y = F.value(x)
    y isa Float64 && isinf(y) && return false
    return true
end

"""
$(TYPEDSIGNATURES)

Return all the equality constraints of `m` as a tuple `(constraint,
bound)` representing `constraint == bound` for each entry.
"""
function equality_constraints(m::C.ConstraintTree)
    sink = Vector{Tuple{Union{LinearValueP,LinearValue},Ex}}()
    get_equality_constraints(m, sink)
    sink
end

"""
$(TYPEDSIGNATURES)

Helper function for [`equality_constraints`](@ref).
"""
function get_equality_constraints(m::C.ConstraintTree, sink)
    get_equality_constraints.(values(m), Ref(sink))
end

"""
$(TYPEDSIGNATURES)

Helper function for [`equality_constraints`](@ref).
"""
function get_equality_constraints(c::C.Constraint, sink)
    if c.bound isa C.EqualToT # TODO check if works
        push!(sink, (C.value(c), Ex(c.bound.equal_to)))
    end
end

"""
$(TYPEDSIGNATURES)

Return all the inequality constraints of `m` as a tuple of bounds converted to
the form `(constraint, upper)` representing `constraint ≤ upper` for each
entry.
"""
function inequality_constraints(m::C.ConstraintTree)
    sink = Vector{Tuple{Union{LinearValueP,LinearValue},Ex,Ex}}()
    get_inequality_constraints(m, sink)

    [
        [ # lower bounds to l ≤ x → -x ≤ -l
            (-val, -lower) for (val, lower, _) in sink if is_constrained(lower)
        ]
        [ # upper bounds are already in the correct format
            (val, upper) for (val, _, upper) in sink if is_constrained(upper)
        ]
    ]
end

"""
$(TYPEDSIGNATURES)

Helper function for [`inequality_constraints`](@ref).
"""
function get_inequality_constraints(m::C.ConstraintTree, sink)
    get_inequality_constraints.(values(m), Ref(sink))
end

"""
$(TYPEDSIGNATURES)

Helper function for [`inequality_constraints`](@ref).
"""
function get_inequality_constraints(c::C.Constraint, sink)
    if c.bound isa C.BetweenT # TODO check if works
        push!(sink, (C.value(c), Ex(c.bound.lower), Ex(c.bound.upper)))
    end
end
