```@meta
EditURL = "2-differentiate-enzyme-model.jl"
```

# Differentiating enzyme constrained metabolic models

````@example 2-differentiate-enzyme-model
using DifferentiableMetabolism
using AbstractFBCModels
using Symbolics
using ConstraintTrees
using COBREXA
using Tulip
using JSONFBCModels
import Downloads: download
using CairoMakie

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

include("../../test/data_static.jl")

model = load_model("e_coli_core.json")
````

unconstrain glucose

````@example 2-differentiate-enzyme-model
rids = [x["id"] for x in model.reactions]
glc_idx = first(indexin(["EX_glc__D_e"], rids))
model.reactions[glc_idx]["lower_bound"] = -1000.0
````

constrain PFL to zero

````@example 2-differentiate-enzyme-model
pfl_idx = first(indexin(["PFL"], rids))
model.reactions[pfl_idx]["upper_bound"] = 0.0

kcats = vcat([Symbolics.@variables $x for x in Symbol.(keys(ecoli_core_reaction_kcats))]...)

rid_kcat = Dict(zip(keys(ecoli_core_reaction_kcats), kcats))

parameter_values =
    Dict(kid => ecoli_core_reaction_kcats[rid] * 3.6 for (rid, kid) in rid_kcat) # change unit to k/h

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

km = build_kinetic_model(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity = capacitylimitation,
)

ec_solution, _, _, _ = optimized_constraints_with_parameters(
    km,
    parameter_values;
    objective = km.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution
````

This solution contains many inactive reactions

````@example 2-differentiate-enzyme-model
sort(abs.(collect(values(ec_solution.fluxes))))
````

And also many inactive gene products.

````@example 2-differentiate-enzyme-model
sort(abs.(collect(values(ec_solution.gene_product_amounts))))
````

With theory, you can show that this introduces flux variability into the
solution, making it non-unique, and consequently non-differentiable. To fix
this, one need to prune the model, to include only the active reactions and
genes. This can be differentiated uniquely. Below we build a pruned kinetic
model, by removing all the reactions, metabolites, and genes that are not
active.

````@example 2-differentiate-enzyme-model
flux_zero_tol = 1e-6 # these bounds make a real difference!
gene_zero_tol = 1e-6

pruned_reaction_isozymes =
    prune_reaction_isozymes(reaction_isozymes, ec_solution, gene_zero_tol)

pruned_gene_product_molar_masses =
    prune_gene_product_molar_masses(gene_product_molar_masses, ec_solution, gene_zero_tol)

pkm = build_kinetic_model(
    prune_model(model, ec_solution, flux_zero_tol, gene_zero_tol);
    reaction_isozymes = pruned_reaction_isozymes,
    gene_product_molar_masses = pruned_gene_product_molar_masses,
    capacity = capacitylimitation,
)

ec_solution2, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    pkm,
    parameter_values;
    objective = pkm.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)
````

Notice, the solution is exactly the same as before, except that all the
inactive elements are gone.

````@example 2-differentiate-enzyme-model
ec_solution2
````

no zero fluxes

````@example 2-differentiate-enzyme-model
sort(abs.(collect(values(ec_solution2.fluxes))))
````

no zero genes

````@example 2-differentiate-enzyme-model
sort(abs.(collect(values(ec_solution2.gene_product_amounts))))
````

prune the kcats, leaving only those that are actually used

````@example 2-differentiate-enzyme-model
pruned_kcats = [
    begin
        x = first(values(v))
        isnothing(x.kcat_forward) ? x.kcat_reverse : x.kcat_forward
    end for v in values(pruned_reaction_isozymes)
]

parameters = [capacitylimitation; kcats]

sens, vids = differentiate(
    pkm,
    pkm.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    parameters;
    scale = true, # unitless sensitivities
)
````

look at glycolysis and oxidative phosphorylation only

````@example 2-differentiate-enzyme-model
subset_ids = [
    r["id"] for r in model.reactions[indexin(string.(keys(pkm.fluxes)), rids)] if
    r["subsystem"] in ["Glycolysis/Gluconeogenesis", "Oxidative Phosphorylation"]
]

flux_idxs = findall(x -> string(last(x)) in subset_ids && first(x) == :fluxes, vids)
flux_ids = last.(vids[flux_idxs])

iso_idxs = findall(x -> string(x[2]) in subset_ids && first(x) != :fluxes, vids)
iso_ids = [v[2] for v in vids[iso_idxs]]

param_idxs = findall(x -> string(x) in subset_ids, parameters)
param_ids = parameters[param_idxs]
````

Flux sensitivities

````@example 2-differentiate-enzyme-model
f, a, hm = heatmap(
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
Colorbar(f[1, 2], hm)
f
````

Isozyme sensitivities. Note, the gene products themselves are not variables in
the formulation of the kinetic model. It inherits its structure from COBREXA,
where the gene products are derived variables. If you want the sensitivities
of the gene products themselves, you just need to multiply the isozyme
sensitivity with the subunit stoichiometry of the relevant gene products.

````@example 2-differentiate-enzyme-model
f, a, hm = heatmap(
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
Colorbar(f[1, 2], hm)
f
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

