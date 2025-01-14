
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
# Construct a community of two interacting E. coli's that share oxoglutarate and
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

include("../../test/data_static.jl")

# Load model, and convert to CanonicalModel for ease of use
model = convert(A.CanonicalModel.Model, X.load_model("e_coli_core.json"))

for rid in ["ACt2r", "ETOHt2r", "PYRt2", "SUCCt3"]
    model.reactions[rid].gene_association_dnf = [["s0001"]]
end

# Modify the model a little bit
model.reactions["EX_glc__D_e"].lower_bound = -1000.0 # capacity bound suffices
model.reactions["PFL"].upper_bound = 0.0 # aerobic simulation

# add a transporter for glutamine
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

# Populate reaction_isozymes with parameters
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

# Add gene product molar mass and capacity constraint info
gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

wt = X.enzyme_constrained_flux_balance_constraints( # reference model, will be used to get some bound info from
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = Dict(
        :total => (
            Symbol.(A.genes(model)),
            C.Between(0, 50),
        ),
    ),
    interface = :identifier_prefixes,
)

km = X.optimized_values( #src
    wt, #src
    optimizer = T.Optimizer, #src
    objective = wt.objective.value, #src
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.8865945023142839, atol = TEST_TOLERANCE) #src

gln_ko = deepcopy(model) # cannot produce glutamine
delete!(gln_ko.reactions, "GLNS") # KO reaction

km = X.enzyme_constrained_flux_balance_analysis( #src
    gln_ko; #src
    reaction_isozymes, #src
    gene_product_molar_masses, #src
    capacity = Dict( #src
        :total => ( #src
            Symbol.(A.genes(model)), #src
            C.Between(0, 50), #src
        ), #src
    ), #src
    optimizer = T.Optimizer, #src
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.0, atol = TEST_TOLERANCE) #src

gln_ko.reactions["EX_gln__L_e"].lower_bound = -1000.0 # open exchange bound to other organism

akg_ko = deepcopy(model) # cannot produce oxoglutarate
delete!(akg_ko.reactions, "ICDHyr") # KO reaction

km = X.enzyme_constrained_flux_balance_analysis( #src
    akg_ko; #src
    reaction_isozymes, #src
    gene_product_molar_masses, #src
    capacity = Dict( #src
        :total => ( #src
            Symbol.(A.genes(model)), #src
            C.Between(0, 50), #src
        ), #src
    ), #src
    optimizer = T.Optimizer, #src
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.0, atol = TEST_TOLERANCE) #src

akg_ko.reactions["EX_akg_e"].lower_bound = -1000.0 # open exchange bound to other organism

ec_gln_ko = X.enzyme_constrained_flux_balance_constraints(
    gln_ko;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = Dict(
        :total => ( 
            Symbol.(A.genes(model)), 
            C.Between(0, 50), 
        ), 
    ), 
    interface = :identifier_prefixes,
)

ec_akg_ko = X.enzyme_constrained_flux_balance_constraints(
    akg_ko;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = Dict(
        :total => ( 
            Symbol.(A.genes(model)), 
            C.Between(0, 50), 
        ), 
    ), 
    interface = :identifier_prefixes,
)

# bound the exchanges - adopt similar bounds to wt wrt to the environment
boundf(id) = begin
    ex_id = first(id)
    wt.interface.exchanges[ex_id].bound
end

# create KO pairs to be joined via their exchanges
ko_pairs = [
    :gln => (ec_gln_ko, ec_gln_ko.interface.exchanges, 0.5) # simulate at abundance of 0.5 for each organism
    :akg => (ec_akg_ko, ec_akg_ko.interface.exchanges, 0.5)
]

# create community model (interface joins them)
x = X.interface_constraints(ko_pairs...; bound = boundf)

# set all growth rates equal
x *= :equalgrowth^C.Constraint(x.gln.objective.value - x.akg.objective.value, C.EqualTo(0))

# set objective as any of the biomass functions (they are constrained equal)
x *= :objective^C.Constraint(x.gln.objective.value, nothing)

sol = X.optimized_values( # test that community can grow
    x;
    optimizer = T.Optimizer,
    objective = x.objective.value,
    sense = X.Maximal,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

@test isapprox(sol.objective, 0.8641022797127501, atol = TEST_TOLERANCE) #src

# now prune the base models
# slighty more involved this time

#md # !!! info "Pruning the model can be tricky"
#md #     Important to ensure no zero fluxes are included in the model - uniqueness.

p_akg_ko, p_r_iso_akg = D.prune_model(
    akg_ko,
    sol.akg.fluxes,
    sol.akg.gene_product_amounts,
    reaction_isozymes,
    sol.akg.isozyme_forward_amounts,
    sol.akg.isozyme_reverse_amounts,
    1e-6,
    1e-6,
);

p_gln_ko, p_r_iso_gln = D.prune_model(
    gln_ko,
    sol.gln.fluxes,
    sol.gln.gene_product_amounts,
    reaction_isozymes,
    sol.gln.isozyme_forward_amounts,
    sol.gln.isozyme_reverse_amounts,
    1e-6,
    1e-6,
);

ec_gln_ko = X.enzyme_constrained_flux_balance_constraints(
    p_gln_ko;
    reaction_isozymes = p_r_iso_gln,
    gene_product_molar_masses,
    capacity = Dict(
        :total => ( 
            Symbol.(A.genes(model)), 
            C.Between(0, 50), 
        ), 
    ), 
    interface = :identifier_prefixes,
)

ec_akg_ko = X.enzyme_constrained_flux_balance_constraints(
    p_akg_ko;
    reaction_isozymes = p_r_iso_akg,
    gene_product_molar_masses,
    capacity = Dict(
        :total => ( 
            Symbol.(A.genes(model)), 
            C.Between(0, 50), 
        ), 
    ), 
    interface = :identifier_prefixes,
)


# bound the exchanges - adopt similar bounds to wt wrt to the environment
env = deepcopy(ec_gln_ko.interface.exchanges)

# delete the exchange of akg and gln to the environment! Will introduce a 0 flux if not
delete!(env, :EX_akg_e)
delete!(env, :EX_gln__L_e)
env

boundf(id) = begin
    ex_id = first(id)
    env[ex_id].bound
end

# ignore the deleted bound, it is still in the interface of each organism and will generate an error if not ignored
ignoref(id, p) = begin
    (id == :gln || id == :akg) && (:EX_gln__L_e in p || :EX_akg_e in p)
end

# investigate the sensitivity of the variables wrt to the abundances of each organism
F.@variables a_gln a_akg

# add parameters into community construction
ko_pairs = [
    :gln => (ec_gln_ko, ec_gln_ko.interface.exchanges, a_gln)
    :akg => (ec_akg_ko, ec_akg_ko.interface.exchanges, a_akg)
]

# create community model (interface joins them)
x = X.interface_constraints(ko_pairs...; bound = boundf, ignore = ignoref)

# join partner exchanges - need to do this to avoid making a 0 env variable
x.interface_balance *=
    :akg^C.Constraint(
        a_akg * x.akg.interface.exchanges.EX_akg_e.value -
        a_gln * x.gln.interface.exchanges.EX_akg_e.value,
        C.EqualTo(0),
    )
x.interface_balance *=
    :gln^C.Constraint(
        a_akg * x.akg.interface.exchanges.EX_gln__L_e.value -
        a_gln * x.gln.interface.exchanges.EX_gln__L_e.value,
        C.EqualTo(0),
    )

# set all growth rates equal
x *= :equalgrowth^C.Constraint(x.gln.objective.value - x.akg.objective.value, C.EqualTo(0))

# set objective as any of the biomass functions (they are constrained equal)
x *= :objective^C.Constraint(x.gln.objective.value, nothing)

param_vals = Dict(:a_gln => 0.5, :a_akg => 0.5)

psol = D.optimized_values(
    x,
    param_vals;
    optimizer = T.Optimizer,
    objective = x.objective.value,
    sense = X.Maximal,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)
psol.tree

@test isapprox(sol.objective, psol.tree.objective, atol = TEST_TOLERANCE) #src

# differentiate wrt to the params
dparams = [:a_akg, :a_gln]

# prepare derivatives
pm_kkt, vids = D.differentiate_prepare_kkt(x, x.objective.value, dparams)

sens = D.differentiate_solution(
    pm_kkt,
    psol.primal_values,
    psol.equality_dual_values,
    psol.inequality_dual_values,
    param_vals,
    scale = true, # unitless sensitivities
)

# only look at how the abundances impact the environmental exchanges
env_exs = string.(last.(vids)[end-7:end])

# lets look at the environmental exchanges
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
