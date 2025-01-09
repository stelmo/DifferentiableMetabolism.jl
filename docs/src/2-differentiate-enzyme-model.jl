
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

# # Differentiating enzyme constrained metabolic models

import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import AbstractFBCModels as A
import JSONFBCModels as JFBC
import COBREXA as X
import Tulip as T
import Clarabel as Q
import Downloads: download
import CairoMakie as CM

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

include("../../test/data_static.jl")

# Load model, and convert to CanonicalModel for ease of use
model = convert(A.CanonicalModel.Model, X.load_model("e_coli_core.json"))

for rid in ["ACt2r", "ETOHt2r", "PYRt2", "SUCCt3"]
    model.reactions[rid].gene_association_dnf = [["s0001"]]
end

# Modify the model a little bit
model.reactions["EX_glc__D_e"].lower_bound = -1000.0 # capacity bound suffices
model.reactions["PFL"].upper_bound = 0.0 # aerobic simulation

# Create parameters of all kcats
rid_kcat = Dict(k => Ex(Symbol(k)) for (k, _) in ecoli_core_reaction_kcats)

# Create a lookup table to map parameters to values
parameter_values = Dict{Symbol,Float64}()

# Create a symbolic reaction_isozyme structure to feed into COBREXA
reaction_isozymes = Dict{String,Dict{String,X.IsozymeT{Ex}}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
float_reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}() #src

# Populate reaction_isozymes with parameters
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available

    for (i, grr) in enumerate(grrs)

        kcat = ecoli_core_reaction_kcats[rid] * 3.6 # change unit to k/h
        parameter_values[Symbol(rid)] = kcat # to substitute later

        d = get!(reaction_isozymes, rid, Dict{String,X.IsozymeT{Ex}}()) # NB: IsozymeT
        d["isozyme_$i"] = X.IsozymeT{Ex}(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = rid_kcat[rid], # assume forward and reverse have the same kcat
            kcat_reverse = rid_kcat[rid],
        )

        d2 = get!(float_reaction_isozymes, rid, Dict{String,X.Isozyme}()) #src
        d2["isozyme_$i"] = X.Isozyme( #src
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), #src
            kcat_forward = kcat, #src
            kcat_reverse = kcat, #src
        ) #src
    end
end

#md # !!! tip "Use the generalized IsozymeT struct from COBREXA"
#md #     Note, COBREXA.jl exports `Isozyme` which is specialized to Float64. To use parameters as shown here, you _must_ use the more general type `IsozymeT`.

# Add gene product molar mass and capacity constraint info
gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

F.@variables capacitylimitation
parameter_values[:capacitylimitation] = 50.0 # mg enzyme/gDW

# Create and solve a COBREXA enzyme constrained model
km = X.enzyme_constrained_flux_balance_constraints( # kinetic model
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = capacitylimitation,
)

ec_solution = D.optimized_constraints_with_parameters(
    km,
    parameter_values;
    objective = km.objective.value,
    optimizer = T.Optimizer,
    modifications = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution.tree

ec_solution_cobrexa = X.enzyme_constrained_flux_balance_analysis( #src
    model; #src
    reaction_isozymes = float_reaction_isozymes, #src
    gene_product_molar_masses = ecoli_core_gene_product_masses, #src
    capacity = 50.0, #src
    optimizer = T.Optimizer, #src
) #src

@test isapprox( #src
    ec_solution.tree.objective, #src
    ec_solution_cobrexa.objective; #src
    atol = TEST_TOLERANCE, #src
) #src

# Note, this solution contains many inactive reactions
sort(collect(ec_solution.tree.fluxes), by = ComposedFunction(abs, last))

@test any(values(ec_solution.tree.fluxes) .≈ 0) #src

# And also many inactive gene products.

sort(collect(ec_solution.tree.gene_product_amounts), by = last)

@test any(isapprox.(values(ec_solution.tree.gene_product_amounts), 0, atol = 1e-8)) #src

# With theory, you can show that this introduces flux variability into the
# solution, making it non-unique, and consequently non-differentiable. To fix
# this, one need to prune the model, to include only the active reactions and
# genes. This can be differentiated uniquely. Below we build a pruned kinetic
# model, by removing all the reactions, metabolites, and genes that are not
# active.

flux_zero_tol = 1e-6 # these bounds make a real difference!
gene_zero_tol = 1e-6

pruned_model, pruned_reaction_isozymes = D.prune_model(
    model,
    ec_solution.tree.fluxes,
    ec_solution.tree.gene_product_amounts,
    reaction_isozymes,
    ec_solution.tree.isozyme_forward_amounts,
    ec_solution.tree.isozyme_reverse_amounts,
    flux_zero_tol,
    gene_zero_tol,
);

pruned_model

@test length(pruned_reaction_isozymes) < length(reaction_isozymes) #src
@test length(pruned_model.reactions) < length(model.reactions) #src
@test length(pruned_model.metabolites) < length(model.metabolites) #src
@test length(pruned_model.genes) < length(model.genes) #src

pkm = X.enzyme_constrained_flux_balance_constraints( # pruned kinetic model
    pruned_model;
    reaction_isozymes = pruned_reaction_isozymes,
    gene_product_molar_masses,
    capacity = [("total", A.genes(pruned_model), capacitylimitation)],
)

pruned_solution = D.optimized_constraints_with_parameters(
    pkm,
    parameter_values;
    objective = pkm.objective.value,
    optimizer = T.Optimizer,
    modifications = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

# Notice, the solution is exactly the same as before, except that all the
# inactive elements are gone.

pruned_solution.tree

# no zero fluxes and all fluxes are made positive!
sort(collect(pruned_solution.tree.fluxes), by = ComposedFunction(abs, last))

# no zero genes
sort(abs.(collect(values(pruned_solution.tree.gene_product_amounts))))

@test isapprox( #src
    pruned_solution.tree.objective, #src
    ec_solution.tree.objective; #src
    atol = TEST_TOLERANCE, #src
) #src

@test all(values(pruned_solution.tree.fluxes) .> 0.0) #src

@test all( #src
    abs(ec_solution.tree.fluxes[k]) - abs(pruned_solution.tree.fluxes[k]) <= 1e-6 for #src
    k in intersect(keys(ec_solution.tree.fluxes), keys(pruned_solution.tree.fluxes)) #src
) #src

# Now we will differentiate the solution. First select parameters that will be differentiated.

kcats = Symbol.(keys(pruned_reaction_isozymes))
parameters = [:capacitylimitation; kcats]

# Next prepare the model for differentiation
pkm_kkt, vids = D.differentiate_prepare_kkt(pkm, pkm.objective.value, parameters)

sens = D.differentiate_solution(
    pkm_kkt,
    pruned_solution.primal_values,
    pruned_solution.equality_dual_values,
    pruned_solution.inequality_dual_values,
    parameter_values,
    scale = true, # unitless sensitivities
)

# look at oxidative phosphorylation only
subset_ids = [:CYTBD, :NADH16, :ATPS4r]

flux_idxs = findall(x -> last(x) in subset_ids && first(x) == :fluxes, vids)
flux_ids = last.(vids[flux_idxs])

iso_idxs = findall(x -> x[2] in subset_ids && occursin("isozyme", string(x[1])), vids)
iso_ids = [v[2] for v in vids[iso_idxs]]

param_idxs = findall(x -> x in subset_ids, parameters)
param_ids = parameters[param_idxs]

# Flux sensitivities
f, a, hm = CM.heatmap(
    sens[flux_idxs, param_idxs]';
    axis = (
        xlabel = "kcat",
        xticks = (1:length(param_ids), string.(param_ids)),
        xticklabelrotation = -pi / 2,
        ylabel = "Flux",
        yticks = (1:length(flux_ids), string.(flux_ids)),
        title = "Flux sensitivities",
    ),
)
CM.Colorbar(f[1, 2], hm)
f

# Isozyme sensitivities. Note, the gene products themselves are not variables in
# the formulation of the kinetic model. It inherits its structure from COBREXA,
# where the gene products are derived variables. If you want the sensitivities
# of the gene products themselves, you just need to multiply the isozyme
# sensitivity with the subunit stoichiometry of the relevant gene products.

f, a, hm = CM.heatmap(
    sens[iso_idxs, param_idxs]';
    axis = (
        xlabel = "kcat",
        xticks = (1:length(param_ids), string.(param_ids)),
        xticklabelrotation = -pi / 2,
        ylabel = "Isozyme",
        yticks = (1:length(iso_ids), string.(iso_ids)),
        title = "Isozyme sensitivities",
    ),
)
CM.Colorbar(f[1, 2], hm)
f

old_atps4r_kcat = parameter_values[:ATPS4r] #src
parameter_values[:ATPS4r] *= 1.0001 #src
kcat_diff = parameter_values[:ATPS4r] - old_atps4r_kcat #src
fin_diff_sol = D.optimized_constraints_with_parameters( #src
    pkm, #src
    parameter_values; #src
    objective = pkm.objective.value, #src
    optimizer = T.Optimizer, #src
    modifications = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
fd_sens = Dict( #src
    k => #src
        (fin_diff_sol.tree.fluxes[k] - pruned_solution.tree.fluxes[k]) / kcat_diff * #src
        old_atps4r_kcat / pruned_solution.tree.fluxes[k] for #src
    k in keys(pruned_solution.tree.fluxes) #src
) #src
fidxs = findall(x -> x[1] == :fluxes, vids) #src
fids = last.(vids)[fidxs] #src
p_atps4r_idx = findfirst(:ATPS4r .== parameters) #src
anal_sens = Dict(y => sens[x, p_atps4r_idx] for (x, y) in zip(fidxs, fids)) #src
@test all(abs(fd_sens[k] - anal_sens[k]) <= TEST_TOLERANCE for k in keys(fd_sens)) #src
