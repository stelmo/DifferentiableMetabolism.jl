
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

# # Parameter estimation using proteomics and flux data

import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import AbstractFBCModels as A
import COBREXA as X
import Tulip as T
import Clarabel as Q
import CairoMakie as CM

# load a small test model
include("../../test/simple_model.jl");

# prune model
delete!(model.reactions, "r5")
delete!(model.genes, "g4")
delete!(model.genes, "g5")
delete!(model.genes, "g3")
model.reactions["r4"].gene_association_dnf = [["g2"]]
model.reactions["r1"].lower_bound = -1000.0
model.reactions["r2"].lower_bound = -1000.0

# now models looks like this

# ![simple_model](./assets/simple_model_pruned.svg)

F.@variables r3 r4

reaction_isozymes = Dict(
    "r3" => Dict(
        "isozyme1" => X.IsozymeT(
            gene_product_stoichiometry = Dict("g1" => 1.0),
            kcat_forward = r3,
            kcat_reverse = nothing,
        ),
    ),
    "r4" => Dict(
        "isozyme1" => X.IsozymeT(
            gene_product_stoichiometry = Dict("g2" => 1.0),
            kcat_forward = r4,
            kcat_reverse = nothing,
        ),
    ),
)

gene_product_molar_masses = Dict("g1" => 20.0, "g2" => 10.0)

F.@variables capacitylimitation

true_parameter_values = Dict(:capacitylimitation => 50.0, :r3 => 2.0, :r4 => 3.0)

km = X.enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = capacitylimitation,
)

sol = D.optimized_constraints_with_parameters(
    km,
    true_parameter_values;
    objective = km.objective.value,
    optimizer = T.Optimizer,
)

sol.tree.fluxes

# create a loss function
measured = [
    sol.tree.fluxes.r1,
    sol.tree.fluxes.r3,
    sol.tree.isozyme_forward_amounts.r3.isozyme1,
    sol.tree.isozyme_forward_amounts.r4.isozyme1,
]

km *=
    :loss^C.Constraint(;
        value = 0.5 * (
            C.squared(km.fluxes.r1.value - measured[1]) +
            C.squared(km.fluxes.r3.value - measured[2]) +
            C.squared(km.isozyme_forward_amounts.r3.isozyme1.value - measured[3]) +
            C.squared(km.isozyme_forward_amounts.r4.isozyme1.value - measured[4])
        ),
        bound = nothing,
    )

estimated_parameters = Dict(:capacitylimitation => 50.0, :r3 => 5.0, :r4 => 1.0) # initial values
η = 0.1 # learning rate

losses = Float64[]

kmKKT, vids =
    D.differentiate_prepare_kkt(km, km.loss.value, [:r3, :r4, :capacitylimitation])

for k = 1:150

    _sol = D.optimized_constraints_with_parameters(
        km,
        estimated_parameters;
        objective = km.loss.value,
        optimizer = Q.Optimizer,
        sense = X.Minimal,
        modifications = [X.silence],
    )
    push!(losses, _sol.tree.loss)

    sens = D.differentiate_solution(
        kmKKT,
        _sol.primal_values,
        _sol.equality_dual_values,
        _sol.inequality_dual_values,
        estimated_parameters,
    )
    measured_idxs = [1, 3, 12, 11]

    x = [
        _sol.tree.fluxes.r1,
        _sol.tree.fluxes.r3,
        _sol.tree.isozyme_forward_amounts.r3.isozyme1,
        _sol.tree.isozyme_forward_amounts.r4.isozyme1,
    ]

    dL_dx = x - measured # derivative of loss function with respect to optimization variables
    dL_dkcats = sens[measured_idxs, :]' * dL_dx # derivative of loss function with respect to parameters

    estimated_parameters[:r3] -= η * dL_dkcats[1]
    estimated_parameters[:r4] -= η * dL_dkcats[2]
end

CM.lines(losses; axis = (xlabel = "Iterations", ylabel = "L2 loss"))


@test abs(estimated_parameters[:r3] - true_parameter_values[:r3]) <= 0.1 #src
@test abs(estimated_parameters[:r4] - true_parameter_values[:r4]) <= 0.1 #src
@test all(losses[2:end] .<= losses[1:end-1]) #src
