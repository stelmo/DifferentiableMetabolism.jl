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
            kcat_backward = rid_kcat[rid],
        )
        break
    end
end

gene_product_molar_masses = Dict(k => v for (k, v) in ecoli_core_gene_product_masses)

Symbolics.@variables capacitylimitation
parameter_values[capacitylimitation] = 50.0 # mg enzyme/gDW

# build parameter model
import ConstraintTrees as C
import COBREXA as X

m = fbc_model_constraints(model)

isozyme_amounts = isozyme_amount_variables(
    Symbol.(keys(reaction_isozymes)),
    rid -> Symbol.(keys(reaction_isozymes[string(rid)])),
)

# allocate variables for everything (nb. += wouldn't associate right here)
m =
    m +
    :fluxes_forward^unsigned_positive_contribution_variables(m.fluxes) +
    :fluxes_backward^unsigned_negative_contribution_variables(m.fluxes) +
    :isozyme_forward_amounts^isozyme_amounts +
    :isozyme_backward_amounts^isozyme_amounts +
    :gene_product_amounts^C.variables(
        keys = Symbol.(AbstractFBCModels.genes(model)),
        bounds = C.Between(0, Inf),
    )

# connect all parts with constraints
m =
    m *
    :directional_flux_balance^sign_split_constraints(
        positive = m.fluxes_forward,
        negative = m.fluxes_backward,
        signed = m.fluxes,
    ) *
    :isozyme_flux_forward_balance^isozyme_flux_constraints(
        m.isozyme_forward_amounts,
        m.fluxes_forward,
        (rid, isozyme) -> X.maybemap(
            x -> x.kcat_forward,
            X.maybeget(reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :isozyme_flux_backward_balance^isozyme_flux_constraints(
        m.isozyme_backward_amounts,
        m.fluxes_backward,
        (rid, isozyme) -> X.maybemap(
            x -> x.kcat_backward,
            X.maybeget(reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :gene_product_isozyme_balance^gene_product_isozyme_constraints(
        m.gene_product_amounts,
        (m.isozyme_forward_amounts, m.isozyme_backward_amounts),
        (rid, isozyme) -> X.maybemap(
            x -> [(Symbol(k), v) for (k, v) in x.gene_product_stoichiometry],
            X.maybeget(reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :gene_product_capacity^C.Constraint(
        value = sum(
            gpa.value * gene_product_molar_masses[String(gp)] for
            (gp, gpa) in m.gene_product_amounts
        ),
        bound = ParameterBetween(0, capacitylimitation),
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
ec_solution.gene_product_amounts

sort(abs.(collect(values(ec_solution.fluxes))))
sort(abs.(collect(values(ec_solution.gene_product_amounts))))

# now prune
zerotol = 1e-7
rids = [string(k) for (k, v) in ec_solution.fluxes if abs(v) > zerotol]
gids = [string(k) for (k, v) in ec_solution.gene_product_amounts if abs(v) > zerotol]
mids = Set(
    mid for rid in rids for
    mid in keys(AbstractFBCModels.reaction_stoichiometry(model, rid))
)

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

# 
m = fbc_model_constraints(pruned)

pruned_reaction_isozymes = Dict(
    k => v for (k, v) in reaction_isozymes if k in AbstractFBCModels.reactions(pruned)
)

isozyme_amounts = isozyme_amount_variables(
    Symbol.(keys(reaction_isozymes)),
    rid -> Symbol.(keys(reaction_isozymes[string(rid)])),
)

# allocate variables for everything (nb. += wouldn't associate right here)
m =
    m +
    :fluxes_forward^unsigned_positive_contribution_variables(m.fluxes) +
    :fluxes_backward^unsigned_negative_contribution_variables(m.fluxes) +
    :isozyme_forward_amounts^isozyme_amounts +
    :isozyme_backward_amounts^isozyme_amounts +
    :gene_product_amounts^C.variables(
        keys = Symbol.(AbstractFBCModels.genes(pruned)),
        bounds = C.Between(0, Inf),
    )

# connect all parts with constraints
m =
    m *
    :directional_flux_balance^sign_split_constraints(
        positive = m.fluxes_forward,
        negative = m.fluxes_backward,
        signed = m.fluxes,
    ) *
    :isozyme_flux_forward_balance^isozyme_flux_constraints(
        m.isozyme_forward_amounts,
        m.fluxes_forward,
        (rid, isozyme) -> X.maybemap(
            x -> x.kcat_forward,
            X.maybeget(pruned_reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :isozyme_flux_backward_balance^isozyme_flux_constraints(
        m.isozyme_backward_amounts,
        m.fluxes_backward,
        (rid, isozyme) -> X.maybemap(
            x -> x.kcat_backward,
            X.maybeget(pruned_reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :gene_product_isozyme_balance^gene_product_isozyme_constraints(
        m.gene_product_amounts,
        (m.isozyme_forward_amounts, m.isozyme_backward_amounts),
        (rid, isozyme) -> X.maybemap(
            x -> [(Symbol(k), v) for (k, v) in x.gene_product_stoichiometry],
            X.maybeget(pruned_reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :gene_product_capacity^C.Constraint(
        value = sum(
            gpa.value * gene_product_molar_masses[String(gp)] for
            (gp, gpa) in m.gene_product_amounts
        ),
        bound = ParameterBetween(0, capacitylimitation),
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

sort(abs.(collect(values(ec_solution.fluxes))))
sort(abs.(collect(values(ec_solution.gene_product_amounts))))

## 
sens = differentiate(
    m,
    m.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    [capacitylimitation; kcats],
    primal_zero_tol = 1e-6,
)

parameters = [capacitylimitation; kcats]
objective = m.objective.value
primal_zero_tol = 1e-8
rank_zero_tol = 1e-8
import LinearAlgebra: :(\), rank, qr
import SparseArrays: sparse, sparsevec, findnz, dropzeros, dropzeros!
