import JSONFBCModels
import Tulip
using Symbolics
using DifferentiableMetabolism
import AbstractFBCModels
using Symbolics
using ConstraintTrees
using COBREXA
using Tulip
using Gurobi

include("data_static.jl")

model = load_model("e_coli_core.json")

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
        d = get!(reaction_isozymes, rid, Dict{String,SimpleIsozyme}())
        d["isozyme_"*string(i)] = ParameterIsozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = rid_kcat[rid],
            kcat_backward = rid_kcat[rid],
        )
        break # single isozyme only
    end
end

Symbolics.@variables capacitylimitation
total_enzyme_capacity = 0.1 * 1000.0 # mg enzyme/gDW
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

ec_solution, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    m,
    parameter_values;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution

ec_solution.fluxes
ec_solution.enzymes

sort(abs.(collect(values(ec_solution.fluxes))))
sort(abs.(collect(values(ec_solution.enzymes))))


# now prune
zerotol = 1e-7
rids = [string(k) for (k, v) in ec_solution.fluxes if abs(v) > zerotol]
gids = [string(k) for (k, v) in ec_solution.enzymes if abs(v) > zerotol]
mids = Set(mid for rid in rids for mid in keys(AbstractFBCModels.reaction_stoichiometry(model, rid)))

d_rids = setdiff(AbstractFBCModels.reactions(model), rids)
d_mids = setdiff(AbstractFBCModels.metabolites(model), mids)
d_gids = setdiff(AbstractFBCModels.genes(model), gids)

pruned = convert(AbstractFBCModels.CanonicalModel.Model, model)

for rid in d_rids
    delete!(pruned.reactions, rid)
end
for mid in d_mids
    delete!(pruned.metabolites, mid)
end
for gid in d_gids
    delete!(pruned.genes, gid)
end

pruned

AbstractFBCModels.n_reactions(pruned)
AbstractFBCModels.n_reactions(model)

AbstractFBCModels.n_metabolites(pruned)
AbstractFBCModels.n_metabolites(model)

AbstractFBCModels.n_genes(pruned)
AbstractFBCModels.n_genes(model)


m = COBREXA.fbc_model_constraints(pruned)
m += :enzymes^COBREXA.enzyme_variables(pruned)
m = COBREXA.add_enzyme_constraints!(m, reaction_isozymes)
m *=
    :total_proteome_bound^enzyme_capacity(
        m.enzymes,
        molar_masses,
        AbstractFBCModels.genes(pruned),
        ParameterBetween(0, capacitylimitation),
    )
m.fluxes.EX_glc__D_e.bound = ConstraintTrees.Between(-1000.0, 0.0) # undo glucose important bound from original model

ec_solution, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    m,
    parameter_values;
    objective = m.objective.value,
    optimizer = Gurobi.Optimizer,
    # modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)
ec_solution.fluxes
ec_solution.enzymes

sort(abs.(collect(values(ec_solution.fluxes))))
sort(abs.(collect(values(ec_solution.enzymes))))

sens = differentiate(
    m,
    m.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    [capacitylimitation; kcats],
    zero_tol = 1e-6,
)

