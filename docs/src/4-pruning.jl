
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import AbstractFBCModels as A
import JSONFBCModels as JFBC
import COBREXA as X
import Tulip as T


include("../../test/data_static.jl")

# Load model, and convert to CanonicalModel for ease of use
model = convert(A.CanonicalModel.Model, X.load_model("e_coli_core.json"))

for rid in ["ACt2r", "ETOHt2r", "PYRt2", "SUCCt3"]
    model.reactions[rid].gene_association_dnf = [["s0001"]]
end

# Modify the model a little bit
model.reactions["EX_glc__D_e"].lower_bound = -1000.0 # capacity bound suffices
model.reactions["PFL"].upper_bound = 0.0 # aerobic simulation

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

gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)


m = X.enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = 50.0,
)

sol = X.optimized_values(
    m,
    objective = m.objective.value,
    optimizer = T.Optimizer,
    settings = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)
