using DifferentiableMetabolism
using COBREXA, ConstraintTrees, FastDifferentiation
using Test

import ConstraintTrees as C

const Expression = FastDifferentiation.Node
const cval = FastDifferentiation.constant_value

t = C.LinearValueT(Expression(0))
@test isempty(t.idxs)
@test isempty(t.weights)

t = C.QuadraticValueT(Expression(0))
@test isempty(t.idxs)
@test isempty(t.weights)

t = C.LinearValueT(Expression(2))
@test all(t.idxs .== [0])
@test all(cval.(t.weights) .== [2])

# L
t = Expression(2) + LinearValue(2)

t = Expression(4) - LinearValue(2)

t = Expression(2) * LinearValue(2)

t = LinearValue(2) / Expression(2)

t = LinearValue(2) + LinearValueT(Expression(2))

t = LinearValue(2) - LinearValueT(Expression(2))

# q

t = Expression(2) + QuadraticValue(2)

t = Expression(4) - QuadraticValue(2)

t = Expression(2) * QuadraticValue(2)

t = QuadraticValue(2) / Expression(2)

t = QuadraticValue(2) + QuadraticValueT(Expression(2))

t = QuadraticValue(2) - QuadraticValueT(Expression(2))

# promotion

t = LinearValueT(Expression(2)) + QuadraticValueT(Expression(2)) 

t = LinearValueT(Expression(2)) - QuadraticValueT(Expression(2)) 

t = LinearValueT(Expression(2)) * LinearValueT(3) 

