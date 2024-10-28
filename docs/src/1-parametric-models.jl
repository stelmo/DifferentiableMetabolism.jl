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

using FastDifferentiation
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
@variables r2bound m3bound

m.fluxes.r2 =
    ConstraintTrees.Constraint(m.fluxes.r2.value, -2 * ParameterBetween(r2bound, 0))

m.flux_stoichiometry.m3 =
    ConstraintTrees.Constraint(m.flux_stoichiometry.m3.value, ParameterEqualTo(m3bound) / 2)

# # add parametric constraints
p = make_variables(:p, 4)

m *=
    :linparam^ConstraintTrees.Constraint(
        value = p[1] * m.fluxes.r1.value + p[2] * m.fluxes.r2.value,
        bound = -ParameterBetween(p[3], 0),
    )

# substitute params into model
parameter_substitutions = Dict(
    :r2bound => 4.0,
    :m3bound => 0.1, # lose some mass here
    :p1 => 1.0,
    :p2 => 1.0,
    :p3 => 4.0,
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
parameter_substitutions[:m3bound] = 0.0

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
)
m_noparams.fluxes

@test isapprox(m_noparams.objective, 4.0; atol = TEST_TOLERANCE) #src

# ## Quadratic parameters also work

q = make_variables(:q, 6)

m.objective = ConstraintTrees.Constraint(
    value = sum(
        rxn.value * rxn.value * qi for (qi, rxn) in zip(collect(q), values(m.fluxes))
    ),
    bound = nothing,
)

m *= :objective_bound^ConstraintTrees.Constraint(value = m.fluxes.r6.value, bound = 2.0)

parameter_substitutions =
    merge(parameter_substitutions, Dict(v.node_value => 1.0 for v in q))

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Clarabel.Optimizer,
    sense = Minimal,
)
m_noparams.fluxes

v(x) = DifferentiableMetabolism.substitute(x, _ -> throw("no subst!")) #src
@test isapprox(m_noparams.objective, 11.0; atol = TEST_TOLERANCE) #src
pb = ParameterBetween(1, 2.0) #src
@test v((-pb).lower) == -1 && v((-pb).upper) == -2 #src
@test v((pb * 2).lower) == 2 && v((pb * 2).upper) == 4 #src
@test v((pb / 2).lower) == 0.5 && v((pb / 2).upper) == 1 #src
pe = ParameterEqualTo(1) #src
@test v((-pe).equal_to) == -1 #src
@test v((pe * 2).equal_to) == 2 #src
@test v((pe / 2).equal_to) == 0.5 #src
plv1 = ParameterLinearValue(ConstraintTrees.LinearValue(1)) #src
@test plv1.idxs == [0] && all(v.(plv1.weights) .== 1.0) #src
plv2 = ParameterLinearValue([0], [1.0]) #src
@test plv2.idxs == [0] && all(v.(plv2.weights) .== 1.0) #src
plv3 = ParameterLinearValue(0) #src
@test isempty(plv3.idxs) && isempty(plv3.weights) #src
plv4 = ParameterLinearValue(1) #src
@test plv4.idxs == [0] && all(v.(plv4.weights) .== 1.0) #src
plv5 = convert(ParameterLinearValue, 1) #src
@test plv5.idxs == [0] && all(v.(plv5.weights) .== 1.0) #src
plv6 = zero(ParameterLinearValue) #src
@test isempty(plv6.idxs) && isempty(plv6.weights) #src
plv7 = plv5 + 1 #src
@test all(v.(plv7.idxs) .== 0) && all(v.(plv7.weights) .== 2.0) #src
plv8 = ParameterLinearValue([1, 2, 3], [1, 2, 3]) #src
plv9 = ParameterLinearValue([1, 3], [4.0, 5.0]) #src
plv10 = plv8 + plv9 #src
@test plv10.idxs == [1, 2, 3] && v.(plv10.weights) == [5, 2, 8] #src
plv10 = plv9 + plv8 #src
@test plv10.idxs == [1, 2, 3] && v.(plv10.weights) == [5, 2, 8] #src
plv11 = plv5 - 2 #src
@test plv11.idxs == [0] && v.(plv11.weights) == [-1] #src
plv12 = 2 - plv5 #src
@test plv12.idxs == [0] && v.(plv12.weights) == [1] #src
plv13 = plv5 / -2 #src
@test plv13.idxs == [0] && v.(plv13.weights) == [-0.5] #src
plv14 = plv8 - plv9 #src
@test plv14.idxs == [1, 2, 3] && v.(plv14.weights) == [-3, 2, -2] #src
plv15 = 1.0 + ParameterLinearValue(0) #src
@test plv15.idxs == [0] && v.(plv15.weights) == [1] #src
pqv1 = ParameterQuadraticValue([(1, 1)], [1]) #src
pqv2 = ParameterQuadraticValue(ConstraintTrees.QuadraticValue([(1, 1)], [1])) #src
@test pqv2.idxs == [(1, 1)] && v.(pqv2.weights) == [1] #src
pqv3 = ParameterQuadraticValue(0) #src
@test isempty(pqv3.idxs) && isempty(pqv3.weights) #src
pqv4 = ParameterQuadraticValue(ParameterLinearValue([1], [1])) #src
@test pqv4.idxs == [(0, 1)] && v.(pqv4.weights) == [1] #src
pqv5 = ParameterQuadraticValue(1) #src
@test pqv5.idxs == [(0, 0)] && v.(pqv5.weights) == [1] #src
pqv6 = convert(ParameterQuadraticValue, 1) #src
@test pqv6.idxs == [(0, 0)] && v.(pqv6.weights) == [1] #src
pqv7 = convert(ParameterQuadraticValue, ParameterLinearValue([1], [1])) #src
@test pqv7.idxs == [(0, 1)] && v.(pqv7.weights) == [1] #src
pqv8 = zero(ParameterQuadraticValue) #src
@test isempty(pqv8.idxs) && isempty(pqv8.weights) #src
pqv9 = 1 + pqv1 #src
@test pqv9.idxs == [(0, 0), (1, 1)] && v.(pqv9.weights) == [1, 1] #src
pqv10 = ParameterLinearValue([1], [1]) + pqv1 #src
@test pqv10.idxs == [(0, 1), (1, 1)] && v.(pqv10.weights) == [1, 1] #src
pqv11 = 3 - pqv1 #src
@test pqv11.idxs == [(0, 0), (1, 1)] && v.(pqv11.weights) == [3, -1] #src
pqv12 = ParameterLinearValue([1], [1]) - pqv1 #src
@test pqv12.idxs == [(0, 1), (1, 1)] && v.(pqv12.weights) == [1, -1] #src
pqv13 = pqv1 - ParameterLinearValue([1], [1]) #src
@test pqv13.idxs == [(0, 1), (1, 1)] && v.(pqv13.weights) == [-1, 1] #src
pqv14 = 42 * pqv1 #src
@test pqv14.idxs == [(1, 1)] && v.(pqv14.weights) == [42] #src
pqv15 = pqv1 - (3 * pqv1) #src
@test pqv15.idxs == [(1, 1)] && v.(pqv15.weights) == [-2] #src
pqv16 = pqv1 / 2 #src
@test pqv16.idxs == [(1, 1)] && v.(pqv16.weights) == [0.5] #src
pqv17 = ParameterLinearValue([1, 2], [2, 1]) * ParameterLinearValue([2], [1]) #src
@test pqv17.idxs == [(1, 2), (2, 2)] && v.(pqv17.weights) == [2, 1] #src
pqv18 = ParameterQuadraticValue([(1, 1), (2, 2)], [1.0, 3.0]) #src
pqv19 = pqv18 + pqv1 #src
@test pqv19.idxs == [(1, 1), (2, 2)] && v.(pqv19.weights) == [2, 3] #src
pqv20 = pqv1 + pqv18 #src
@test pqv20.idxs == [(1, 1), (2, 2)] && v.(pqv19.weights) == [2, 3] #src
pqv21 = ParameterQuadraticValue([(1, 1), (2, 2)], [1, 2]) #src
@test pqv21.idxs == [(1, 1), (2, 2)] && v.(pqv21.weights) == [1, 2] #src
pqv22 = convert(ParameterQuadraticValue, ConstraintTrees.LinearValue([1], [1])) #src
@test pqv22.idxs == [(0, 1)] && v.(pqv22.weights) == [1] #src
pqv23 = convert( #src
    ParameterQuadraticValue, #src
    ConstraintTrees.QuadraticValue([(1, 1), (2, 2)], [1, 2]), #src
) #src
@test pqv23.idxs == [(1, 1), (2, 2)] && v.(pqv23.weights) == [1, 2] #src
pqv24 = 5 - ConstraintTrees.QuadraticValue([(1, 1)], [1]) #src
@test pqv24.idxs == [(0, 0), (1, 1)] && v.(pqv24.weights) == [5, -1] #src
pqv25 = pqv18 - 5.0 #src
@test pqv25.idxs == [(0, 0), (1, 1), (2, 2)] && v.(pqv25.weights) == [-5, 1, 3] #src
pr1 = 1 + ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test pr1.idxs == [0, 1, 2] && v.(pr1.weights) == [1, 1, 2] #src
pr2 = ConstraintTrees.LinearValue([1, 2], [1, 2]) - 3 #src
@test pr2.idxs == [0, 1, 2] && v.(pr2.weights) == [-3, 1, 2] #src
pr2b = 3 - ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test pr2b.idxs == [0, 1, 2] && v.(pr2b.weights) == [3, -1, -2] #src
pr3 = 3 * ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test pr3.idxs == [1, 2] && v.(pr3.weights) == [3, 6] #src
pr3 = ConstraintTrees.LinearValue([1, 2], [1, 2]) / 2 #src
@test pr3.idxs == [1, 2] && v.(pr3.weights) == [0.5, 1] #src
pr4 = ConstraintTrees.LinearValue([1, 2], [1, 2]) + ParameterLinearValue([1, 2], [1, 2]) #src
@test pr4.idxs == [1, 2] && v.(pr4.weights) == [2, 4] #src
pr5 = ParameterLinearValue([1, 2], [2, 3]) - ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test pr5.idxs == [1, 2] && v.(pr5.weights) == [1, 1] #src
pr6 = 1 + ConstraintTrees.QuadraticValue([(1, 1)], [1]) #src
@test pr6.idxs == [(0, 0), (1, 1)] && v.(pr6.weights) == [1, 1] #src
pr7 = ConstraintTrees.QuadraticValue([(1, 1)], [1]) - 1 #src
@test pr7.idxs == [(0, 0), (1, 1)] && v.(pr7.weights) == [-1, 1] #src
pr8 = 2 * ConstraintTrees.QuadraticValue([(1, 1)], [1]) #src
@test pr8.idxs == [(1, 1)] && v.(pr8.weights) == [2] #src
pr9 = ConstraintTrees.QuadraticValue([(1, 1)], [1]) / 2 #src
@test pr9.idxs == [(1, 1)] && v.(pr9.weights) == [0.5] #src
pr10 = #src
    ConstraintTrees.QuadraticValue([(1, 1)], [1]) + ParameterQuadraticValue([(1, 1)], [1]) #src
@test pr10.idxs == [(1, 1)] && v.(pr10.weights) == [2] #src
pr11 = #src
    ParameterQuadraticValue([(1, 1)], [2]) - ConstraintTrees.QuadraticValue([(1, 1)], [1]) #src
@test pr11.idxs == [(1, 1)] && v.(pr11.weights) == [1] #src
pr12 = ConstraintTrees.LinearValue([1, 2], [1, 2]) + ParameterQuadraticValue([(1, 1)], [2]) #src
@test pr12.idxs == [(0, 1), (1, 1), (0, 2)] && v.(pr12.weights) == [1, 2, 2] #src
pr13 = ParameterQuadraticValue([(1, 1)], [2]) - ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test pr13.idxs == [(0, 1), (1, 1), (0, 2)] && v.(pr13.weights) == [-1, 2, -2] #src
pr14 = ParameterLinearValue([1, 2], [2, 3]) * ConstraintTrees.LinearValue([1, 2], [1, 2]) #src
@test pr14.idxs == [(1, 1), (1, 2), (2, 2)] && v.(pr14.weights) == [2, 7, 6] #src
pr15 = ConstraintTrees.LinearValue([1, 2], [1, 2]) - ParameterQuadraticValue([(1, 1)], [2]) #src
@test pr15.idxs == [(0, 1), (1, 1), (0, 2)] && v.(pr15.weights) == [1, -2, 2] #src
@test v( #src
    ConstraintTrees.substitute( #src
        ParameterQuadraticValue([(0, 0), (1, 1)], [1, 2]), #src
        FastDifferentiation.Node.([1.0, 2.0]), #src
    ), #src
) == 3 #src
