import JSONFBCModels
import Tulip
using Symbolics
using DifferentiableMetabolism
import AbstractFBCModels
using Symbolics
using ConstraintTrees
using COBREXA
using Tulip

include("data_static.jl")

model = load_model("e_coli_core.json")

Symbolics.@variables kcats[1:length(ecoli_core_reaction_kcats)]
rid_kcat = Dict(zip(keys(ecoli_core_reaction_kcats), kcats))
parameter_values = Dict(kid => ecoli_core_reaction_kcats[rid] * 3.6/100.0 for (rid, kid) in rid_kcat) # k/h

reaction_isozymes = Dict{String,Dict{String,ParameterIsozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for rid in AbstractFBCModels.reactions(model)
    grrs = AbstractFBCModels.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, rid, Dict{String,SimpleIsozyme}())
        d["isozyme_"*string(i)] = ParameterIsozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = rid_kcat[rid],
            kcat_backward = rid_kcat[rid],
        )
    end
end

Symbolics.@variables capacitylimitation
total_enzyme_capacity = 0.1 * 1000.0 * 100 # mg enzyme/gDW
parameter_values[capacitylimitation] = total_enzyme_capacity

molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

# Build differentiable enzyme constrained model
m = COBREXA.fbc_model_constraints(model)
m += :enzymes^COBREXA.enzyme_variables(model)
m = COBREXA.add_enzyme_constraints!(m, reaction_isozymes)
m *=
    :total_proteome_bound^enzyme_capacity(
        m.enzymes,
        molar_masses,
        AbstractFBCModels.genes(model),
        ParameterBetween(0, capacitylimitation),
    )
m.fluxes.EX_glc__D_e.bound = ConstraintTrees.Between(-1000.0, 0.0) # undo glucose important bound from original model

# add small quadratic weight
m.objective = ConstraintTrees.Constraint(
    value =  1e-4*(sum(x.value * x.value for x in values(m.fluxes)) + sum(x.value * x.value for x in values(m.enzymes))) - m.objective.value,
    bound = nothing,
) 

using Gurobi
ec_solution, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    m,
    parameter_values;
    objective = m.objective.value,
    optimizer = Gurobi.Optimizer,
    sense = COBREXA.Minimal,
)

ec_solution.fluxes.BIOMASS_Ecoli_core_w_GAM
ec_solution.fluxes
ec_solution.enzymes

sens = differentiate(
    m,
    m.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    [capacitylimitation; kcats],
)

objective = m.objective.value
parameters = [capacitylimitation; kcats]
