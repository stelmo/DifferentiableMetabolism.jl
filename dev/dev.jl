using COBREXA, JSON, HiGHS, FastDifferentiation, DifferentiableMetabolism
import AbstractFBCModels.CanonicalModel as CM
import AbstractFBCModels as A
import JSONFBCModels


mass_bound = 350.0
optimizer = HiGHS.Optimizer
gene_zero_tol = 1e-8
flux_zero_tol = 1e-8


model = convert(CM.Model, load_model("dev/Alloascoidea_hylecoeti.json"));
reaction_isozymes = Dict(
    k => Dict(
        "$(k)_$i" => Isozyme(
            d["gene_product_stoichiometry"],
            d["kcat_forward"],
            d["kcat_reverse"],
        ) for (i, d) in enumerate(v)
    ) for (k, v) in JSON.parsefile("dev/Alloascoidea_hylecoeti_isozymes.json")
);
gene_product_molar_masses =
    Dict(x => y for (x, y) in JSON.parsefile("dev/Alloascoidea_hylecoeti_molar_mass.json"))
capacity = [("uncategorized", A.genes(model), mass_bound)]
# open exchange reactions
for (r, rxn) in model.reactions
    if contains(rxn.name, "exchange")
        model.reactions[r].lower_bound = rxn.lower_bound * 1000.0
        model.reactions[r].upper_bound = rxn.upper_bound * 1000.0
    end
end

ec_solution = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer,
)

reaction_parameter_isozymes = Dict(
    r => Dict(
        k => ParameterIsozyme(
            v.gene_product_stoichiometry,
            FastDifferentiation.Node(Symbol("$(k)_f")),
            FastDifferentiation.Node(Symbol("$(k)_r")),
        ) for (k, v) in isos
    ) for (r, isos) in reaction_isozymes
);

pruned_reaction_isozymes =
    prune_reaction_isozymes(reaction_parameter_isozymes, ec_solution, gene_zero_tol);

pruned_gene_product_molar_masses =
    prune_gene_product_molar_masses(gene_product_molar_masses, ec_solution, gene_zero_tol);

pruned_model = prune_model(model, ec_solution, flux_zero_tol, gene_zero_tol)
#make reactions uni-directional
for (r,rxn) in pruned_model.reactions
    if ec_solution.fluxes[r] < 0 
        pruned_model.reactions[r].upper_bound = 0.0
    else
        pruned_model.reactions[r].lower_bound = 0.0
    end
end


pkm = build_kinetic_model(
    pruned_model;
    reaction_isozymes=pruned_reaction_isozymes,
    gene_product_molar_masses=pruned_gene_product_molar_masses,
    capacity,
);

parameter_values = Dict{Symbol,Float64}();
for (r, iso) in pruned_reaction_isozymes
    for (k, v) in iso
        parameter_values[Symbol(k, :_f)] = reaction_isozymes[r][k].kcat_forward
        parameter_values[Symbol(k, :_r)] = reaction_isozymes[r][k].kcat_reverse
    end
end
ec_solution2, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    pkm,
    parameter_values;
    objective=pkm.objective.value,
    optimizer=HiGHS.Optimizer,
    modifications=[silence],
);
parameters = collect(keys(parameter_values));

sens, vids = differentiate(
    pkm,
    pkm.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    parameters;
    scale=true, # unitless sensitivities
);
