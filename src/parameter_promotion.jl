
#=
Copyright (c) 2023, Heinrich-Heine University Duesseldorf
Copyright (c) 2023, University of Luxembourg

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

# Promote LinearValue to ParameterLinearValue if it interacts with a parameter

# Nums and LinearValues
Base.:+(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = b + a
Base.:+(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = a + ParameterLinearValue(b)

Base.:-(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = -b + a
Base.:-(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = a - ParameterLinearValue(b)

Base.:*(a::Symbolics.Num, b::ConstraintTrees.LinearValue) = b * a
Base.:*(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = ParameterLinearValue(a.idxs, b .* a.weights)

Base.:/(a::ConstraintTrees.LinearValue, b::Symbolics.Num) = ParameterLinearValue(a.idxs, a.weights ./ b)

# LinearValue and ParameterLinearValues
Base.:+(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) = b + a
Base.:+(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) = a + ParameterLinearValue(b)

Base.:-(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) = -b + a
Base.:-(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) = a - ParameterLinearValue(b)

# Promote QuadraticValue to ParameterQuadraticValue if it interacts with a parameter

# Nums and QuadraticValues
Base.:+(a::Symbolics.Num, b::ConstraintTrees.QuadraticValue) = b + a
Base.:+(a::ConstraintTrees.QuadraticValue, b::Symbolics.Num) = a + ParameterQuadraticValue(b)

Base.:-(a::Symbolics.Num, b::ConstraintTrees.QuadraticValue) = -b + a
Base.:-(a::ConstraintTrees.QuadraticValue, b::Symbolics.Num) = a - ParameterQuadraticValue(b)

Base.:*(a::Symbolics.Num, b::ConstraintTrees.QuadraticValue) = b * a
Base.:*(a::ConstraintTrees.QuadraticValue, b::Symbolics.Num) = ParameterQuadraticValue(a.idxs, b .* a.weights)

Base.:/(a::ConstraintTrees.QuadraticValue, b::Symbolics.Num) = ParameterQuadraticValue(a.idxs, a.weights ./ b)

# QuadraticValues and ParameterQuadraticValues
Base.:+(a::ConstraintTrees.QuadraticValue, b::ParameterQuadraticValue) = b + a
Base.:+(a::ParameterQuadraticValue, b::ConstraintTrees.QuadraticValue) = a + ParameterQuadraticValue(b)

Base.:-(a::ConstraintTrees.QuadraticValue, b::ParameterQuadraticValue) = -b + a
Base.:-(a::ParameterQuadraticValue, b::ConstraintTrees.QuadraticValue) = a - ParameterQuadraticValue(b)

# Cross interaction terms

# LinearValues and ParameterQuadraticValues (+, - only)
Base.:+(a::ConstraintTrees.LinearValue, b::ParameterQuadraticValue) = b + a
Base.:+(a::ParameterQuadraticValue, b::ConstraintTrees.LinearValue) = a + ParameterQuadraticValue(b)

Base.:-(a::ConstraintTrees.LinearValue, b::ParameterQuadraticValue) = -b + a
Base.:-(a::ParameterQuadraticValue, b::ConstraintTrees.LinearValue) = a - ParameterQuadraticValue(b)

# LinearValue * ParameterLinearValue
Base.:*(a::ParameterLinearValue, b::ConstraintTrees.LinearValue) = b * a
Base.:*(a::ConstraintTrees.LinearValue, b::ParameterLinearValue) = ParameterLinearValue(a) * b
