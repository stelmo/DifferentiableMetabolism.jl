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
const Expression = FastDifferentiation.Node
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

m.fluxes.r2 = ConstraintTrees.Constraint(m.fluxes.r2.value, BetweenT(-2 * r2bound, Expression(0)))

m.flux_stoichiometry.m3 = ConstraintTrees.Constraint(m.flux_stoichiometry.m3.value, EqualToT(m3bound) / 2)

# # add parametric constraints
p = make_variables(:p, 4)

m *= :linparam^ConstraintTrees.Constraint(
        value = p[1] * m.fluxes.r1.value + p[2] * m.fluxes.r2.value,
        bound = BetweenT(-p[3], Expression(0)),
    )

# substitute params into model
parameter_substitutions = Dict(
    :r2bound => 4.0,
    :m3bound => 0.1, # lose some mass here
    :p1 => 1.0,
    :p2 => 1.0,
    :p3 => 4.0,
)

m_substituted = DifferentiableMetabolism.substitute(m, k -> parameter_substitutions[k])

optimized_values(m_substituted,objective = m.objective.value, optimizer = Tulip.Optimizer,)

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m_substituted,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
)

# alternatively, solve directly
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

