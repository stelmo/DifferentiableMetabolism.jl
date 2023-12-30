# Copyright (c) 2023, Heinrich-Heine University Duesseldorf                #src
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

# # Differentiating enzyme constrained metabolic models 

using DifferentiableMetabolism
import AbstractFBCModels
using Symbolics
using ConstraintTrees
using COBREXA
using Tulip

include(joinpath("..", "test", "simple_model.jl")) #hide

# ![simple_model](./assets/simple_model.svg)

model

# ## Add enzyme kinetic information

Symbolics.@variables kcats_forward[1:4] kcats_backward[1:4]

reaction_isozymes = Dict(
    "r3" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[1], kcats_backward[1]),
    ),
    "r4" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[2], kcats_backward[2]),
        "iso2" =>
            ParameterIsozyme(Dict("g2" => 1), kcats_forward[3], kcats_backward[3]),
    ),
    "r5" => Dict(
        "iso1" => ParameterIsozyme(
            Dict("g3" => 1, "g4" => 2),
            kcats_forward[4],
            kcats_backward[4],
        ),
    ),
)

gene_molar_masses = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0, "g5" => 1.0)

# ## Add a parameterized capacity limitation

Symbolics.@variables capacitylimitation

# Build differentiable enzyme constrained model
m = COBREXA.fbc_model_constraints(model)
m += :enzymes^COBREXA.enzyme_variables(model)
m = COBREXA.add_enzyme_constraints!(m, reaction_isozymes)
m *=
    :total_proteome_bound^ConstraintTrees.Constraint(
        value = sum(
            m.enzymes[Symbol(gid)].value * gene_molar_masses[gid] for
            gid in AbstractFBCModels.genes(model)
        ),
        bound = ParameterBetween(0.0, capacitylimitation),
    )
m.enzymes.g3.bound = ConstraintTrees.Between(0, 0.01)
m.fluxes.r2.bound = ConstraintTrees.Between(0, 100)

# substitute params into model
parameters = Dict(
    kcats_forward[1] => 1.0,
    kcats_forward[2] => 2.0,
    kcats_forward[3] => 3.0,
    kcats_forward[4] => 70.0,
    kcats_backward[1] => 1.0,
    kcats_backward[2] => 2.0,
    kcats_backward[3] => 3.0,
    kcats_backward[4] => 70.0,
    capacitylimitation => 0.5,
)

ec_solution = optimized_constraints_with_parameters(
    m,
    parameters;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

@test isapprox(ec_solution.objective, 3.181818181753438, atol = TEST_TOLERANCE)
@test isapprox(ec_solution.enzymes.g4, 0.09090909090607537, atol = TEST_TOLERANCE)
@test isapprox(ec_solution.:total_proteome_bound, 0.5, atol = TEST_TOLERANCE)

ec_solution.fluxes
ec_solution.enzymes

# ## Prune model

reaction_isozymes = Dict(
    "r3" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[1], kcats_backward[1]),
    ),
    "r4" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[2], kcats_backward[2]),
    ),
    "r5" => Dict(
        "iso1" => ParameterIsozyme(
            Dict("g3" => 1, "g4" => 2),
            kcats_forward[4],
            kcats_backward[4],
        ),
    ),
)

delete!(model.genes, "g2")
delete!(model.genes, "g5")

# Build differentiable enzyme constrained model
m = COBREXA.fbc_model_constraints(model)
m += :enzymes^COBREXA.enzyme_variables(model)
m = COBREXA.add_enzyme_constraints!(m, reaction_isozymes)
m *=
    :total_proteome_bound^ConstraintTrees.Constraint(
        value = sum(
            m.enzymes[Symbol(gid)].value * gene_molar_masses[gid] for
            gid in AbstractFBCModels.genes(model)
        ),
        bound = ParameterBetween(0.0, capacitylimitation),
    )
m.enzymes.g3.bound = ConstraintTrees.Between(0, 0.01)
m.fluxes.r2.bound = ConstraintTrees.Between(0, 100)

# substitute params into model
parameters = Dict(
    kcats_forward[1] => 1.0,
    kcats_forward[2] => 2.0,
    kcats_forward[3] => 3.0,
    kcats_forward[4] => 70.0,
    kcats_backward[1] => 1.0,
    kcats_backward[2] => 2.0,
    kcats_backward[3] => 3.0,
    kcats_backward[4] => 70.0,
    capacitylimitation => 0.5,
)

ec_solution = optimized_constraints_with_parameters(
    m,
    parameters;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution.fluxes
ec_solution.enzymes

objective = m.objective

_x, _eq_duals, _ineq_duals = optimized_constraints_with_parameters(
    m,
    parameters;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
    primal_duals_only = true,
)
