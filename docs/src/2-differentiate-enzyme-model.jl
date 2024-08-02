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
using AbstractFBCModels
using Symbolics
using ConstraintTrees
using COBREXA
using Tulip
using JSONFBCModels
using Gurobi
import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

include("./test/data_static.jl")

model = load_model("e_coli_core.json")
# unconstrain glucose
rids = [x["id"] for x in model.reactions]
glc_idx = first(indexin(["EX_glc__D_e"], rids))
model.reactions[glc_idx]["lower_bound"] = -1000.0

Symbolics.@variables kcats[1:length(ecoli_core_reaction_kcats)]
rid_kcat = Dict(zip(keys(ecoli_core_reaction_kcats), kcats))
parameter_values =
    Dict(kid => ecoli_core_reaction_kcats[rid] * 3.6 for (rid, kid) in rid_kcat) # k/h

reaction_isozymes = Dict{String,Dict{String,ParameterIsozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for rid in AbstractFBCModels.reactions(model)
    grrs = AbstractFBCModels.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, rid, Dict{String,ParameterIsozyme}())
        d["isozyme_"*string(i)] = ParameterIsozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = rid_kcat[rid],
            kcat_reverse = rid_kcat[rid],
        )
    end
end

gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

Symbolics.@variables capacitylimitation
parameter_values[capacitylimitation] = 50.0 # mg enzyme/gDW

m = build_kinetic_model(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = capacitylimitation,
)

ec_solution, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_values;
    objective = m.objective.value,
    optimizer = Gurobi.Optimizer,
    # modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution

@test isapprox(ec_solution.objective, 0.7069933828497013; atol = TEST_TOLERANCE)

sort(abs.(collect(values(ec_solution.fluxes)))) # lots of zeros

sort(abs.(collect(values(ec_solution.gene_product_amounts))))

flux_zero_tol = 1e-8
gene_zero_tol = 1e-8

pkm = build_kinetic_model(
    prune_model(model, ec_solution, flux_zero_tol, gene_zero_tol);
    reaction_isozymes = prune_reaction_isozymes(
        reaction_isozymes,
        ec_solution,
        flux_zero_tol,
        gene_zero_tol,
    ),
    gene_product_molar_masses = prune_gene_product_molar_masses(
        gene_product_molar_masses,
        ec_solution,
        flux_zero_tol,
        gene_zero_tol,
    ),
    capacity = capacitylimitation,
)

ec_solution2, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    pkm,
    parameter_values;
    objective = pkm.objective.value,
    optimizer = Gurobi.Optimizer,
    # modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution2

@test isapprox(ec_solution2.objective, 0.7069933828497013; atol = TEST_TOLERANCE)

# no zero fluxes
sort(abs.(collect(values(ec_solution2.fluxes))))

# no zero genes
sort(abs.(collect(values(ec_solution2.gene_product_amounts))))

## 
parameters = [capacitylimitation; kcats]

function tp()
sens, vids = differentiate(
    pkm,
    pkm.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    parameters;
    scale = true,
)
end

@b tp()

