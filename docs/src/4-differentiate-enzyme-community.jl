
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

# # Differentiating enzyme constrained community models
# Construct a community of two interacting E. coli cells that share oxoglutarate and
# glutamine.

import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import AbstractFBCModels as A
import JSONFBCModels as JFBC
import COBREXA as X
import Tulip as T
import Downloads: download
import CairoMakie as CM

# ## Constructing and solving an enzyme constrained community model

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

# Add a transporter for glutamine
model.reactions["GLNt"] = A.CanonicalModel.Reaction(;
    name = "Glutamine transporter",
    lower_bound = -1000.0,
    upper_bound = 1000.0,
    stoichiometry = Dict(
        "gln__L_e" => -1.0,
        "gln__L_c" => 1.0,
        "h_e" => -1.0,
        "h_c" => 1.0,
    ),
)

reaction_isozymes = Dict{String,Dict{String,X.Isozyme}}()

# Populate `reaction_isozymes` with parameters
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available

    for (i, grr) in enumerate(grrs)

        kcat = ecoli_core_reaction_kcats[rid] * 3.6 # change unit to k/h

        d = get!(reaction_isozymes, rid, Dict{String,X.Isozyme}())
        d["isozyme_$i"] = X.Isozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))),
            kcat_forward = kcat,
            kcat_reverse = kcat,
        )
    end
end

# Make the wildtype model `wt`, where we add gene product molar mass and capacity
# constraint info
gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

wt = X.enzyme_constrained_flux_balance_constraints( # reference model, will be used to get some bound info from
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 50.0,
    interface = :identifier_prefixes,
)

km = X.optimized_values( #src
    wt, #src
    optimizer = T.Optimizer, #src
    objective = wt.objective.value, #src
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.8865945023142839, atol = TEST_TOLERANCE) #src

# Make a glutamine auxotrophic mutant `gln_ko`
gln_ko = deepcopy(model) # cannot produce glutamine
delete!(gln_ko.reactions, "GLNS") # KO reaction

km = X.enzyme_constrained_flux_balance_analysis( #src
    gln_ko; #src
    reaction_isozymes, #src
    gene_product_molar_masses, #src
    capacity = 50.0, #src
    optimizer = T.Optimizer, #src
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.0, atol = TEST_TOLERANCE) #src

gln_ko.reactions["EX_gln__L_e"].lower_bound = -1000.0 # open exchange bound to other organism

# Make an oxoglutarate auxotrophic mutant `akg_ko`
akg_ko = deepcopy(model) # cannot produce oxoglutarate
delete!(akg_ko.reactions, "ICDHyr") # KO reaction

km = X.enzyme_constrained_flux_balance_analysis( #src
    akg_ko; #src
    reaction_isozymes, #src
    gene_product_molar_masses, #src
    capacity = 50.0, #src
    optimizer = T.Optimizer, #src
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.0, atol = TEST_TOLERANCE) #src

akg_ko.reactions["EX_akg_e"].lower_bound = -1000.0 # open exchange bound to other organism

# Add isozymes, gene product molar masses, and capacity constraints to the KO mutants
ec_gln_ko = X.enzyme_constrained_flux_balance_constraints(
    gln_ko;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 50.0,
    interface = :identifier_prefixes,
)

ec_akg_ko = X.enzyme_constrained_flux_balance_constraints(
    akg_ko;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 50.0,
    interface = :identifier_prefixes,
)

# Bound the exchanges - adopt similar bounds to wt wrt to the environment
boundf(id) = begin
    ex_id = first(id)
    wt.interface.exchanges[ex_id].bound
end

# Create KO pairs to be joined via their exchanges
ko_pairs = [
    :gln => (ec_gln_ko, ec_gln_ko.interface.exchanges, 0.5) # simulate at abundance of 0.5 for each organism
    :akg => (ec_akg_ko, ec_akg_ko.interface.exchanges, 0.5)
]

# Create community model, where `interface` joins the two enzyme constrained KO models
x = X.interface_constraints(ko_pairs...; bound = boundf)

# Set all growth rates equal
x *= :equalgrowth^C.Constraint(x.gln.objective.value - x.akg.objective.value, C.EqualTo(0))

# Set objective as any of the biomass functions (they are constrained equal)
x *= :objective^C.Constraint(x.gln.objective.value, nothing)

sol = X.optimized_values( # test that community can grow
    x;
    optimizer = T.Optimizer,
    objective = x.objective.value,
    sense = X.Maximal,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

@test isapprox(sol.objective, 0.8641022797127501, atol = TEST_TOLERANCE) #src

# ## Pruning the community

# Now prune the base models. Here it is slighty more involved than in the case of a
# single enzyme constrained model.

#md # !!! info "Pruning the model can be tricky"
#md #     For a unique, and therefore differentiable solution, it is important to ensure that no zero fluxes are
#md #     included in the model.

pruned_akg_ko, pruned_isozymes_akg = D.prune_model(
    akg_ko,
    sol.akg.fluxes,
    sol.akg.gene_product_amounts,
    reaction_isozymes,
    sol.akg.isozyme_forward_amounts,
    sol.akg.isozyme_reverse_amounts,
    1e-6,
    1e-6,
);

pruned_gln_ko, pruned_isozymes_gln = D.prune_model(
    gln_ko,
    sol.gln.fluxes,
    sol.gln.gene_product_amounts,
    reaction_isozymes,
    sol.gln.isozyme_forward_amounts,
    sol.gln.isozyme_reverse_amounts,
    1e-6,
    1e-6,
);

# Make the pruned enzyme constrained KO models
ec_gln_ko = X.enzyme_constrained_flux_balance_constraints(
    pruned_gln_ko;
    reaction_isozymes = pruned_isozymes_gln,
    gene_product_molar_masses,
    capacity = 50.0,
    interface = :identifier_prefixes,
)

ec_akg_ko = X.enzyme_constrained_flux_balance_constraints(
    pruned_akg_ko;
    reaction_isozymes = pruned_isozymes_akg,
    gene_product_molar_masses,
    capacity = 50.0,
    interface = :identifier_prefixes,
)


# Bound the exchanges - adopt similar bounds to wt wrt to the environment
env = deepcopy(ec_gln_ko.interface.exchanges)

#md # !!! info "Delete the KO exchanges"
#md #     Exchanges of akg and gln to the environment need to be deleted! Otherwise, a zero flux will be introduced

delete!(env, :EX_akg_e)
delete!(env, :EX_gln__L_e)
env

boundf(id) = begin
    ex_id = first(id)
    env[ex_id].bound
end

# Ignore the deleted bound, it is still in the interface of each organism and will generate an error if not ignored
ignoref(id, p) = begin
    (id == :gln || id == :akg) && (:EX_gln__L_e in p || :EX_akg_e in p)
end

# ## Investigating sensitivity of variables to the organism abundances
# We want to see how sensitive the variables of both KO models are to the abundances of
# each organism in the community.
F.@variables abundance_gln abundance_akg

# Add parameters into a community construction
ko_pairs = [
    :gln => (ec_gln_ko, ec_gln_ko.interface.exchanges, abundance_gln)
    :akg => (ec_akg_ko, ec_akg_ko.interface.exchanges, abundance_akg)
]

# Create community model (`interface`` joins them)
x = X.interface_constraints(ko_pairs...; bound = boundf, ignore = ignoref)

# Join partner exchanges - need to do this to avoid making a 0 env variable
x.interface_balance *=
    :akg^C.Constraint(
        abundance_akg * x.akg.interface.exchanges.EX_akg_e.value -
        abundance_gln * x.gln.interface.exchanges.EX_akg_e.value,
        C.EqualTo(0),
    )
x.interface_balance *=
    :gln^C.Constraint(
        abundance_akg * x.akg.interface.exchanges.EX_gln__L_e.value -
        abundance_gln * x.gln.interface.exchanges.EX_gln__L_e.value,
        C.EqualTo(0),
    )

# Set all growth rates equal
x *= :equalgrowth^C.Constraint(x.gln.objective.value - x.akg.objective.value, C.EqualTo(0))

# Set objective as any of the biomass functions (they are constrained equal)
x *= :objective^C.Constraint(x.gln.objective.value, nothing)

param_vals = Dict(:abundance_gln => 0.5, :abundance_akg => 0.5)

pruned_sol = D.optimized_values(
    x,
    param_vals;
    optimizer = T.Optimizer,
    objective = x.objective.value,
    sense = X.Maximal,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)
pruned_sol.tree

@test isapprox(sol.objective, pruned_sol.tree.objective, atol = TEST_TOLERANCE) #src

# Want to differenciate with respect to the parameters `abundance_akg` and `abundance_gln`, the respective
# abundances of the two KO mutants
dparams = [:abundance_akg, :abundance_gln]

# Prepare and calculate derivatives
pm_kkt, vids = D.differentiate_prepare_kkt(x, x.objective.value, dparams)

sens = D.differentiate_solution(
    pm_kkt,
    pruned_sol.primal_values,
    pruned_sol.equality_dual_values,
    pruned_sol.inequality_dual_values,
    param_vals,
    scale = true, # unitless sensitivities
)

# Only look at how the abundances impact the environmental exchanges
env_exs = string.(last.(vids)[end-7:end])

# Now we plot the environmental exchanges
fig, ax, hm = CM.heatmap(
    sens[end-7:end, :];
    axis = (
        xticks = (1:8, env_exs),
        xticklabelrotation = -pi / 2,
        yticks = (1:2, string.(dparams)),
    ),
)
CM.Colorbar(fig[1, 2], hm)
fig
