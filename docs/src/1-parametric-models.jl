
# Copyright (c) 2025, Heinrich-Heine University Duesseldorf                #src
# Copyright (c) 2025, University of Luxembourg                             #src
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

import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import COBREXA as X
import Tulip as T
import Clarabel as Q

# ## Load and solve a simple model

# Load a small test model
include("../../test/simple_model.jl");

# ![simple_model](./assets/simple_model.svg)

model

# Build a basic `ConstraintTree` model without parameters
m = X.flux_balance_constraints(model)

# Solve normally
base_model = X.optimized_values(m; optimizer = T.Optimizer, objective = m.objective.value)
base_model.fluxes

@test isapprox(base_model.objective, 1.0; atol = TEST_TOLERANCE) #src

# ## Add parameters to the model

# Make bound of `r2` and mass balance of `m3` parameters
F.@variables r2bound m3bound

m.fluxes.r2 = C.Constraint(m.fluxes.r2.value, C.BetweenT(-2 * r2bound, Ex(0)))

m.flux_stoichiometry.m3 =
    C.Constraint(m.flux_stoichiometry.m3.value, C.EqualToT(m3bound) / 2)

#md # !!! tip "Use the generalized bounds from ConstraintTrees"
#md #     Note, ConstraintTrees.jl exports `Between` and `EqualTo` which are specialized to Float64. To use parameters as shown here, you _must_ use the more general types `BetweenT` and `EqualToT`. Appropriate overloads have been added to simplify type promotion when adding floaty bounds to symbolic bounds.

# # Add parametric constraints
p = F.make_variables(:p, 4)

m *=
    :linparam^C.Constraint(
        value = p[1] * m.fluxes.r1.value + p[2] * m.fluxes.r2.value,
        bound = C.BetweenT(-p[3], Ex(0)),
    )

#md # !!! tip "Use the generalized value types from ConstraintTrees"
#md #     Note, ConstraintTrees.jl exports `LinearValue` and `QuadraticValue` which are specialized to Float64. To use parameters as shown here, you _must_ use the more general types `LinearValueT` and `QuadraticValueT`. Appropriate overloads have been added to simplify type construction and promotion (as used above). But note that `m.linparam.value` is a `ConstraintTrees.LinearValueT{FastDifferentiation.Node}`.

# Substitute parameters into model to yield a "normal" constraint tree model
parameter_substitutions = Dict(
    :r2bound => 4.0,
    :m3bound => 0.1, # lose some mass here
    :p1 => 1.0,
    :p2 => 1.0,
    :p3 => 4.0,
)

m_substituted = D.substitute(m, k -> parameter_substitutions[k])

# This can be solved like any constraint tree
m_normal = X.optimized_values(
    m_substituted,
    objective = m.objective.value,
    optimizer = T.Optimizer,
)

# Alternatively, a convenience function can take care of the substitutions for you
m_noparams = D.optimized_values(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = T.Optimizer,
)

@test isapprox(m_noparams.tree.objective, 3.899999999938411; atol = TEST_TOLERANCE) #src

# ## Change the parameters and re-solve

# Substitute parameters into model
parameter_substitutions[:m3bound] = 0.0

m_noparams2 = D.optimized_values(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = T.Optimizer,
)

m_noparams2.tree.fluxes

@test isapprox(m_noparams2.tree.objective, 4.0; atol = TEST_TOLERANCE) #src

# ## Quadratic parameters also work

q = F.make_variables(:q, 6)

m.objective = C.Constraint(
    value = sum(
        rxn.value * rxn.value * qi for (qi, rxn) in zip(collect(q), values(m.fluxes))
    ),
    bound = nothing,
)

m *= :objective_bound^C.Constraint(value = m.fluxes.r6.value, bound = 2.0)

parameter_substitutions =
    merge(parameter_substitutions, Dict(v.node_value => 1.0 for v in q))

m_noparams3 = D.optimized_values(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Q.Optimizer,
    sense = X.Minimal,
)

m_noparams3.tree.fluxes

@test isapprox(m_noparams3.tree.objective, 11; atol = TEST_TOLERANCE) #src
