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
# parameter_values = Dict(kid => ecoli_core_reaction_kcats[rid]*3600.0 for (rid, kid) in rid_kcat)

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
        break
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

#src these values should be unique (glucose transporter is the only way to get carbon into the system)
@test isapprox(ec_solution.objective, 1.671357282901553, atol = TEST_TOLERANCE) #src
@test isapprox(ec_solution.total_proteome_bound, 0.1, atol = TEST_TOLERANCE) #src
@test isapprox(ec_solution.fluxes.EX_glc__D_e, -49.92966287110028, atol = 0.1) #src
@test isapprox(ec_solution.enzymes.b2417, 0.11859224858442563, atol = 1e-7) #src

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
    primal_zero_tol = 1e-6,
)

