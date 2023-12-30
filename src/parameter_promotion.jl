
# Promote LinearValue to ParameterLinearValue if it interacts with a parameter
Base.:+(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = ParameterLinearValue(a) + b

Base.:+(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = a + ParameterLinearValue(b)

Base.:-(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = ParameterLinearValue(a) - b

Base.:-(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = a - ParameterLinearValue(b)

Base.:*(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = b * a

Base.:*(a::ConstraintTrees.LinearValue, b::Symbolics.Num) =
    ParameterLinearValue(a.idxs, b .* a.weights)

Base.:/(a::ConstraintTrees.LinearValue, b::Symbolics.Num) =
    ParameterLinearValue(a.idxs, a.weights ./ b)

Base.:+(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) =
    ParameterLinearValue(a.idxs, a.weights) + b

Base.:+(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) =
    a + ParameterLinearValue(b.idxs, b.weights)

Base.:-(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) =
    ParameterLinearValue(a.idxs, a.weights) - b
    
Base.:-(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) =
    a - ParameterLinearValue(b.idxs, b.weights)

# Promote QuadraticValue to ParameterQuadraticValue if it interacts with a parameter
