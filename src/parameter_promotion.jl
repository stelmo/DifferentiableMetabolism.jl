
# Promote LinearValue to ParameterLinearValue if it interacts with a parameter

# Nums and LinearValues
Base.:+(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = ParameterLinearValue(a) + b

Base.:+(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = a + ParameterLinearValue(b)

Base.:-(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = ParameterLinearValue(a) - b

Base.:-(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = a - ParameterLinearValue(b)

Base.:*(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = b * a

Base.:*(a::ConstraintTrees.LinearValue, b::Symbolics.Num) =
    ParameterLinearValue(a.idxs, b .* a.weights)

Base.:/(a::ConstraintTrees.LinearValue, b::Symbolics.Num) =
    ParameterLinearValue(a.idxs, a.weights ./ b)

# LinearValue and ParameterLinearValues
Base.:+(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) =
    ParameterLinearValue(a) + b

Base.:+(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) =
    a + ParameterLinearValue(b)

Base.:-(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) =
    ParameterLinearValue(a) - b

Base.:-(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) =
    a - ParameterLinearValue(b)

# Promote QuadraticValue to ParameterQuadraticValue if it interacts with a parameter

# Nums and QuadraticValues
Base.:+(a::Symbolics.Num, b::ConstraintTrees.QuadraticValue) =
    ParameterQuadraticValue(a) + b

Base.:+(a::ConstraintTrees.QuadraticValue, b::Symbolics.Num) =
    a + ParameterQuadraticValue(b)

Base.:-(a::Symbolics.Num, b::ConstraintTrees.QuadraticValue) =
    ParameterQuadraticValue(a) - b

Base.:-(a::ConstraintTrees.QuadraticValue, b::Symbolics.Num) =
    a - ParameterQuadraticValue(b)

Base.:*(a::Symbolics.Num, b::ConstraintTrees.QuadraticValue) = b * a

Base.:*(a::ConstraintTrees.QuadraticValue, b::Symbolics.Num) =
    ParameterQuadraticValue(a.idxs, b .* a.weights)

Base.:/(a::ConstraintTrees.QuadraticValue, b::Symbolics.Num) =
    ParameterQuadraticValue(a.idxs, a.weights ./ b)

# QuadraticValues and ParameterQuadraticValues
Base.:+(a::ConstraintTrees.QuadraticValue, b::ParameterQuadraticValue) =
    ParameterQuadraticValue(a) + b

Base.:+(a::ParameterQuadraticValue, b::ConstraintTrees.QuadraticValue) =
    a + ParameterQuadraticValue(b)

Base.:-(a::ConstraintTrees.QuadraticValue, b::ParameterQuadraticValue) =
    ParameterQuadraticValue(a) - b

Base.:-(a::ParameterQuadraticValue, b::ConstraintTrees.QuadraticValue) =
    a - ParameterQuadraticValue(b)

# Cross interaction terms

# LinearValues and ParameterQuadraticValues (+, - only)
Base.:+(a::ConstraintTrees.LinearValue, b::ParameterQuadraticValue) =
    ParameterQuadraticValue(a) + b

Base.:+(a::ParameterQuadraticValue, b::ConstraintTrees.LinearValue) =
    a + ParameterQuadraticValue(b)

Base.:-(a::ConstraintTrees.LinearValue, b::ParameterQuadraticValue) =
    ParameterQuadraticValue(a) - b

Base.:-(a::ParameterQuadraticValue, b::ConstraintTrees.LinearValue) =
    a - ParameterQuadraticValue(b)

# LinearValue * ParameterLinearValue
Base.:*(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) =
    a * ParameterLinearValue(b)

Base.:*(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) =
    ParameterLinearValue(a) * b
