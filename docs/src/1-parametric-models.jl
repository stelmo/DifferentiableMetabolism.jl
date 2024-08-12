# Copyright (c) 2024, Heinrich-Heine University Duesseldorf                #src
#                                                                          #src
# Licensed under the Apache License, Version 2.0 (the "License");          #src
# you may not use this file except in compliance with the License.         #src
# You may obtain a copy of the License at                                  #src
#                                                                          #src
#     http://www.apache.org/licenses/LICENSE-2.0                           #src
#                                                                          #src
# Unless required by applicable law or agreed to in writing, software      #src
# distributed under the License is distributed on an "AS IS" BASIS,        #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #src
# See the License for the specific language governing permissions and      #src
# limitations under the License.                                           #src

# # Parametric constraint-based metabolic models

using DifferentiableMetabolism

using Symbolics
using ConstraintTrees
using COBREXA
using Tulip
using Clarabel

# ## Load and solve a simple model

# The code used to construct the model is located in `test/simple_model.jl`, but
# it is not shown here for brevity. Below is a visualization of the model.

include("../../test/simple_model.jl"); #hide

# ![simple_model](./assets/simple_model.svg)

model

# Build a basic ConstraintTree model without parameters
m = COBREXA.flux_balance_constraints(model)

# Solve normally
base_model =
    COBREXA.optimized_values(m; optimizer = Tulip.Optimizer, objective = m.objective.value)
base_model.fluxes

@test isapprox(base_model.objective, 1.0; atol = TEST_TOLERANCE) #src

# ## Add parameters to the model

# Make bound of r2 and mass balance of m3 parameters
Symbolics.@variables r2bound m3bound

m.fluxes.r2 =
    ConstraintTrees.Constraint(m.fluxes.r2.value, -2 * ParameterBetween(r2bound, 0))

m.flux_stoichiometry.m3 =
    ConstraintTrees.Constraint(m.flux_stoichiometry.m3.value, ParameterEqualTo(m3bound) / 2)

# # add parametric constraints
Symbolics.@variables p[1:4]

m *=
    :linparam^ConstraintTrees.Constraint(
        value = p[1] * m.fluxes.r1.value + p[2] * m.fluxes.r2.value,
        bound = -ParameterBetween(p[3], 0),
    )

# substitute params into model
parameter_substitutions = Dict(
    r2bound => 4.0,
    m3bound => 0.1, # lose some mass here
    p[1] => 1.0,
    p[2] => 1.0,
    p[3] => 4.0,
)

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
)
m_noparams.fluxes

@test isapprox(m_noparams.objective, 3.899999999938411; atol = TEST_TOLERANCE) #src

# ## Change the parameters and re-solve

# substitute params into model
parameter_substitutions[m3bound] = 0.0

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
)
m_noparams.fluxes

@test isapprox(m_noparams.objective, 4.0; atol = TEST_TOLERANCE) #src

# ## Quadratic parameters also work

Symbolics.@variables q[1:6]

m.objective = ConstraintTrees.Constraint(
    value = sum(
        rxn.value * rxn.value * qi for (qi, rxn) in zip(collect(q), values(m.fluxes))
    ),
    bound = nothing,
)

m *= :objective_bound^ConstraintTrees.Constraint(value = m.fluxes.r6.value, bound = 2.0)

parameter_substitutions = merge(parameter_substitutions, Dict(zip(q, fill(1.0, 6))))

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Clarabel.Optimizer,
    sense = Minimal,
)
m_noparams.fluxes

@test isapprox(m_noparams.objective, 11.0; atol = TEST_TOLERANCE) #src

#src miscellaneous tests
pb = ParameterBetween(1, 2.0) #src
@test (-pb).lower == -1 && (-pb).upper == -2  #src
@test (pb * 2).lower == 2 && (pb * 2).upper == 4  #src
@test (pb / 2).lower == 0.5 && (pb / 2).upper == 1  #src
pe = ParameterEqualTo(1) #src
@test (-pe).equal_to == -1 #src
@test (pe * 2).equal_to == 2 #src
@test (pe / 2).equal_to == 0.5 #src

plv1 = ParameterLinearValue(ConstraintTrees.LinearValue(1)) #src
@test all(plv1.idxs .== [0]) && all(plv1.weights .== [Num(1.0)]) #src
plv2 = ParameterLinearValue([0], [1.0]) #src
@test all(plv2.idxs .== [0]) && all(plv2.weights .== [Num(1.0)]) #src
plv3 = ParameterLinearValue(0) #src
@test isempty(plv3.idxs) && isempty(plv3.weights) #src
plv4 = ParameterLinearValue(1) #src
@test all(plv4.idxs .== [0]) && all(plv4.weights .== [Num(1.0)]) #src
plv5 = convert(ParameterLinearValue, 1) #src
@test all(plv5.idxs .== [0]) && all(plv5.weights .== [Num(1.0)]) #src
plv6 = zero(ParameterLinearValue) #src
@test isempty(plv6.idxs) && isempty(plv6.weights) #src
plv7 = plv5 + 1 #src
@test all(plv7.idxs .== [0]) && all(plv7.weights .== [Num(2.0)]) #src
plv8 = ParameterLinearValue([1, 2, 3], [1, 2, 3]) #src
plv9 = ParameterLinearValue([1, 3], [4.0, 5.0]) #src
plv10 = plv8 + plv9 #src
@test all(plv10.idxs .== [1, 2, 3]) && all(plv10.weights .== [Num(5), Num(2), Num(8)]) #src
plv10 = plv9 + plv8 #src
@test all(plv10.idxs .== [1, 2, 3]) && all(plv10.weights .== [Num(5), Num(2), Num(8)]) #src
plv11 = plv5 - 2 #src
@test all(plv11.idxs .== [0]) && all(plv11.weights .== [Num(-1)]) #src
plv12 = 2 - plv5 #src
@test all(plv12.idxs .== [0]) && all(plv12.weights .== [Num(1)]) #src
plv13 = plv5 / -2 #src
@test all(plv13.idxs .== [0]) && all(plv13.weights .== [Num(-0.5)]) #src
plv14 = plv8 - plv9 #src
@test all(plv14.idxs .== [1, 2, 3]) && all(plv14.weights .== [Num(-3), Num(2), Num(-2)]) #src
plv15 = 1.0 + ParameterLinearValue(0) #src
@test all(plv15.idxs .== [0]) && all(plv15.weights .== [Num(1)]) #src

pqv1 = ParameterQuadraticValue([(1, 1)], [1]) #src
pqv2 = ParameterQuadraticValue(ConstraintTrees.QuadraticValue([(1, 1)], [1])) #src
@test all(pqv2.idxs[1] .== (1, 1)) && all(pqv2.weights .== Num(1)) #src
pqv3 = ParameterQuadraticValue(0) #src
@test isempty(pqv3.idxs) && isempty(pqv3.weights) #src
pqv4 = ParameterQuadraticValue(ParameterLinearValue([1], [1])) #src
@test all(pqv4.idxs .== [(0, 1)]) && all(pqv4.weights .== [Num(1)]) #src
pqv5 = ParameterQuadraticValue(1) #src
@test all(pqv5.idxs .== [(0, 0)]) && all(pqv5.weights .== [Num(1)]) #src
pqv6 = convert(ParameterQuadraticValue, 1) #src
@test all(pqv6.idxs .== [(0, 0)]) && all(pqv6.weights .== [Num(1)]) #src
pqv7 = convert(ParameterQuadraticValue, ParameterLinearValue([1], [1])) #src
@test all(pqv7.idxs .== [(0, 1)]) && all(pqv7.weights .== [Num(1)]) #src
pqv8 = zero(ParameterQuadraticValue) #src
@test isempty(pqv8.idxs) && isempty(pqv8.weights) #src
pqv9 = 1 + pqv1 #src
@test all(pqv9.idxs .== [(0, 0), (1, 1)]) && all(pqv9.weights .== [Num(1), Num(1)]) #src
pqv10 = ParameterLinearValue([1], [1]) + pqv1 #src
@test all(pqv10.idxs .== [(0, 1), (1, 1)]) && all(pqv10.weights .== [Num(1), Num(1)]) #src
pqv11 = 3 - pqv1 #src
@test all(pqv11.idxs .== [(0, 0), (1, 1)]) && all(pqv11.weights .== [Num(3), Num(-1)]) #src
pqv12 = ParameterLinearValue([1], [1]) - pqv1 #src
@test all(pqv12.idxs .== [(0, 1), (1, 1)]) && all(pqv12.weights .== [Num(1), Num(-1)]) #src
pqv13 = pqv1 - ParameterLinearValue([1], [1]) #src
@test all(pqv13.idxs .== [(0, 1), (1, 1)]) && all(pqv13.weights .== [Num(-1), Num(1)]) #src
pqv14 = 42 * pqv1 #src
@test all(pqv14.idxs .== [(1, 1)]) && all(pqv14.weights .== [Num(42)]) #src
pqv15 = pqv1 - (3 * pqv1) #src
@test all(pqv15.idxs .== [(1, 1)]) && all(pqv15.weights .== [Num(-2)]) #src
pqv16 = pqv1 / 2 #src
@test all(pqv16.idxs .== [(1, 1)]) && all(pqv16.weights .== [Num(0.5)]) #src
pqv17 = ParameterLinearValue([1, 2], [2, 1]) * ParameterLinearValue([2], [1]) #src 
@test all(pqv17.idxs .== [(1, 2), (2, 2)]) && all(pqv17.weights .== [Num(2), Num(1)]) #src
pqv18 = ParameterQuadraticValue([(1, 1), (2, 2)], [1.0, 3.0]) #src
pqv19 = pqv18 + pqv1 #src
@test all(pqv19.idxs .== [(1, 1), (2, 2)]) && all(pqv19.weights .== [Num(2), Num(3)]) #src
pqv20 = pqv1 + pqv18 #src
@test all(pqv20.idxs .== [(1, 1), (2, 2)]) && all(pqv19.weights .== [Num(2), Num(3)]) #src
pqv21 = ParameterQuadraticValue([(1, 1), (2, 2)], [1, 2]) #src
@test all(pqv21.idxs .== [(1, 1), (2, 2)]) && all(pqv21.weights .== [Num(1), Num(2)]) #src
pqv22 = convert(ParameterQuadraticValue, ConstraintTrees.LinearValue([1], [1])) #src
@test all(pqv22.idxs .== [(0, 1)]) && all(pqv22.weights .== [Num(1)]) #src
pqv23 = convert( #src
    ParameterQuadraticValue, #src
    ConstraintTrees.QuadraticValue([(1, 1), (2, 2)], [1, 2]), #src
) #src
@test all(pqv23.idxs .== [(1, 1), (2, 2)]) && all(pqv23.weights .== [Num(1), Num(2)]) #src
pqv24 = Num(5) - ConstraintTrees.QuadraticValue([(1, 1)], [1]) #src
@test all(pqv24.idxs .== [(0, 0), (1, 1)]) && all(pqv24.weights .== [Num(5), Num(-1)]) #src
pqv25 = pqv18 - 5.0 #src
@test all(pqv25.idxs .== [(0, 0), (1, 1), (2, 2)]) && #src
      all(pqv25.weights .== [Num(-5), Num(1), Num(3)]) #src
pr1 = Num(1) + ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test all(pr1.idxs .== [0, 1, 2]) && all(pr1.weights .== [Num(1), Num(1), Num(2)]) #src
pr2 = ConstraintTrees.LinearValue([1, 2], [1, 2]) - Num(3) #src
@test all(pr2.idxs .== [0, 1, 2]) && all(pr2.weights .== [Num(-3), Num(1), Num(2)]) #src
pr2b = Num(3) - ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test all(pr2b.idxs .== [0, 1, 2]) && all(pr2b.weights .== [Num(3), Num(-1), Num(-2)]) #src
pr3 = Num(3) * ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test all(pr3.idxs .== [1, 2]) && all(pr3.weights .== [Num(3), Num(6)]) #src
pr3 = ConstraintTrees.LinearValue([1, 2], [1, 2]) / Num(2) #src
@test all(pr3.idxs .== [1, 2]) && all(pr3.weights .== [Num(0.5), Num(1)]) #src
pr4 = ConstraintTrees.LinearValue([1, 2], [1, 2]) + ParameterLinearValue([1, 2], [1, 2]) #src
@test all(pr4.idxs .== [1, 2]) && all(pr4.weights .== [Num(2), Num(4)]) #src
pr5 = ParameterLinearValue([1, 2], [2, 3]) - ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test all(pr5.idxs .== [1, 2]) && all(pr5.weights .== [Num(1), Num(1)]) #src
pr6 = Num(1) + ConstraintTrees.QuadraticValue([(1, 1)], [1]) #src
@test all(pr6.idxs .== [(0, 0), (1, 1)]) && all(pr6.weights .== [Num(1), Num(1)]) #src
pr7 = ConstraintTrees.QuadraticValue([(1, 1)], [1]) - Num(1) #src
@test all(pr7.idxs .== [(0, 0), (1, 1)]) && all(pr7.weights .== [Num(-1), Num(1)]) #src
pr8 = Num(2) * ConstraintTrees.QuadraticValue([(1, 1)], [1]) #src
@test all(pr8.idxs .== [(1, 1)]) && all(pr8.weights .== [Num(2)]) #src
pr9 = ConstraintTrees.QuadraticValue([(1, 1)], [1]) / Num(2)  #src
@test all(pr9.idxs .== [(1, 1)]) && all(pr9.weights .== [Num(0.5)]) #src
pr10 = #src
    ConstraintTrees.QuadraticValue([(1, 1)], [1]) + ParameterQuadraticValue([(1, 1)], [1]) #src
@test all(pr10.idxs .== [(1, 1)]) && all(pr10.weights .== [Num(2)]) #src
pr11 = #src
    ParameterQuadraticValue([(1, 1)], [2]) - ConstraintTrees.QuadraticValue([(1, 1)], [1]) #src
@test all(pr11.idxs .== [(1, 1)]) && all(pr11.weights .== [Num(1)]) #src
pr12 = ConstraintTrees.LinearValue([1, 2], [1, 2]) + ParameterQuadraticValue([(1, 1)], [2]) #src
@test all(pr12.idxs .== [(0, 1), (1, 1), (0, 2)]) && #src
      all(pr12.weights .== [Num(1), Num(2), Num(2)]) #src
pr13 = ParameterQuadraticValue([(1, 1)], [2]) - ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test all(pr13.idxs .== [(0, 1), (1, 1), (0, 2)]) && #src
      all(pr13.weights .== [Num(-1), Num(2), Num(-2)]) #src
pr14 = ParameterLinearValue([1, 2], [2, 3]) * ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test all(pr14.idxs .== [(1, 1), (1, 2), (2, 2)]) && #src
      all(pr14.weights .== [Num(2), Num(7), Num(6)]) #src
pr15 = ConstraintTrees.LinearValue([1, 2], [1, 2]) - ParameterQuadraticValue([(1, 1)], [2]) #src
@test all(pr15.idxs .== [(0, 1), (1, 1), (0, 2)]) && #src
      all(pr15.weights .== [Num(1), Num(-2), Num(2)]) #src
@test ConstraintTrees.substitute( #src
    ParameterQuadraticValue([(0, 0), (1, 1)], [1, 2]), #src
    [Num(1), Num(2)], #src
) == 3 #src
