
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

model.reactions["GLNt"] = A.CanonicalModel.Reaction(;
    name = "Glutamine transporter",
    lower_bound = -1000.0,
    upper_bound = 1000.0,
    stoichiometry = Dict(
        "gln__L_e" => -1.0,
        "gln__L_c" => 1.0,
        "h_e" => -1.0,
        "h_c" => 1.0,
    )
)


reaction_isozymes = Dict{String,Dict{String, X.Isozyme}}() 

# Populate reaction_isozymes with parameters
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available

    for (i, grr) in enumerate(grrs)

        kcat = ecoli_core_reaction_kcats[rid] * 3.6 # change unit to k/h

        d = get!(reaction_isozymes, rid, Dict{String, X.Isozyme}())
        d["isozyme_$i"] = X.Isozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))),
            kcat_forward = kcat, 
            kcat_reverse = kcat, 
        ) 
    end
end

# Add gene product molar mass and capacity constraint info
gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

wt = X.enzyme_constrained_flux_balance_constraints( # reference model
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 50.0,
    interface = :identifier_prefixes
)

km = X.optimized_values( #src
    wt,
    optimizer = T.Optimizer, #src
    objective = wt.objective.value,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.8865945023142839, atol=TEST_TOLERANCE) #src

glutamine_ko = deepcopy(model) # cannot produce glutamine
delete!(glutamine_ko.reactions, "GLNS") # KO reaction

km = X.enzyme_constrained_flux_balance_analysis( #src
    glutamine_ko; #src
    reaction_isozymes, #src
    gene_product_molar_masses, #src
    capacity = 50.0, #src
    optimizer = T.Optimizer, #src
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.0, atol=TEST_TOLERANCE) #src

glutamine_ko.reactions["EX_gln__L_e"].lower_bound = -1000.0 # open exchange bound to other organism


oxoglutarate_ko = deepcopy(model) # cannot produce oxoglutarate
delete!(oxoglutarate_ko.reactions, "ICDHyr") # KO reaction

km = X.enzyme_constrained_flux_balance_analysis( #src
    oxoglutarate_ko; #src
    reaction_isozymes, #src
    gene_product_molar_masses, #src
    capacity = 50.0, #src
    optimizer = T.Optimizer, #src
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
) #src
@test isapprox(km.objective, 0.0, atol=TEST_TOLERANCE) #src

oxoglutarate_ko.reactions["EX_akg_e"].lower_bound = -1000.0 # open exchange bound to other organism

ec_gln_ko = X.enzyme_constrained_flux_balance_constraints(
    glutamine_ko;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 50.0,
    interface = :identifier_prefixes
)


ec_akg_ko = X.enzyme_constrained_flux_balance_constraints(
    oxoglutarate_ko;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 50.0,
    interface = :identifier_prefixes
)

# bound the exchanges - adopt similar bounds to wt wrt to the environment
boundf(id) = begin
    ex_id = first(id)
    wt.interface.exchanges[ex_id].bound
end

# create KO pairs to be joined via their exchanges
ko_pairs = [
    :gln => (ec_gln_ko, ec_gln_ko.interface.exchanges, 0.5)
    :akg => (ec_akg_ko, ec_akg_ko.interface.exchanges, 0.5)
]

# create community model (interface joins them)
x = X.interface_constraints(ko_pairs...; bound = boundf)

# set all growth rates equal
x *= :equalgrowth^C.Constraint(
    x.gln.objective.value - x.akg.objective.value,
    C.EqualTo(0)
)

# set objective as any of the biomass functions (they are constrained equal)
x *= :objective^C.Constraint(
    x.gln.objective.value,
    nothing,
)

sol = X.optimized_values(
    x;
    optimizer = T.Optimizer,
    objective = x.objective.value,
    sense = X.Maximal,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)], #src
)


