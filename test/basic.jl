
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

# make testing easier for values
function value_test(v, idxs, weights)
    @test all(v.idxs .== idxs)
    @test all(F.value.(v.weights) .== weights)
end

# construction

@test isempty(C.LinearValueT(Ex(0)).idxs)
@test isempty(C.LinearValueT(Ex(0)).weights)
@test isempty(C.QuadraticValueT(Ex(0)).idxs)
@test isempty(C.QuadraticValueT(Ex(0)).weights)

value_test(C.LinearValueT(Ex(1)), [0], [1])
value_test(C.QuadraticValueT(Ex(1)), [(0, 0)], [1])

# Linear
t = Ex(2) + C.LinearValue([2], [3])
value_test(t, [0, 2], [2, 3])

t = Ex(4) - C.LinearValue([1], [2])
value_test(t, [0, 1], [4, -2])

t = Ex(2) * C.LinearValue(2)
value_test(t, [0], [4])

t = C.LinearValue([2], [2]) / Ex(2)
value_test(t, [2], [1])

t = C.LinearValue([1], [2]) + C.LinearValueT([2], [Ex(2)])
value_test(t, [1, 2], [2, 2])

t = C.LinearValue([1], [2]) - C.LinearValueT([2], [Ex(2)])
value_test(t, [1, 2], [2, -2])

# Quadratic

t = Ex(2) + C.QuadraticValue(2)
value_test(t, [(0, 0)], [4])

t = Ex(4) - C.QuadraticValue([(0, 1)], [2])
value_test(t, [(0, 0), (0, 1)], [4, -2.0])

t = Ex(2) * C.QuadraticValue(2)
value_test(t, [(0, 0)], [4])

t = C.QuadraticValue(2) / Ex(2)
value_test(t, [(0, 0)], [1])

t = C.QuadraticValue([(0, 1)], [2]) + C.QuadraticValueT([(0, 1)], [Ex(2)])
value_test(t, [(0, 1)], [4.0])

t = C.QuadraticValue(2) - C.QuadraticValueT(Ex(5))
value_test(t, [(0, 0)], [-3.0])

# interaction terms

t = C.LinearValueT(Ex(2)) + C.QuadraticValueT(Ex(3))
value_test(t, [(0, 0)], [5.0])

t = C.LinearValueT(Ex(2)) - C.QuadraticValueT(Ex(3))
value_test(t, [(0, 0)], [-1.0])

t = C.LinearValueT(Ex(2)) * C.LinearValue(3)
value_test(t, [(0, 0)], [6.0])
