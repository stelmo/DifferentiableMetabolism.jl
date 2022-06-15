# DifferentiableMetabolism

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stelmo.github.io/DifferentiableMetabolism.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stelmo.github.io/DifferentiableMetabolism.jl/dev)
[![Build Status](https://github.com/stelmo/DifferentiableMetabolism.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/stelmo/DifferentiableMetabolism.jl/actions/workflows/CI.yml?query=branch%3Amaster)

# Differentiating a simple model

``` 
using DifferentiableMetabolism
using CPLEX
using COBREXA

# Create model
model = StandardModel("SmallModel")
m1 = Metabolite("m1")
m2 = Metabolite("m2")
m3 = Metabolite("m3")
m4 = Metabolite("m4")
m5 = Metabolite("m5")
m6 = Metabolite("m6")

@add_reactions! model begin
    "r1", nothing → m1, 0, 100
    "r2", nothing → m2, 0, 100
    "r3", m1 + m2 → m3, 0, 100
    "r4", m3 → m4 + m5, 0, 100
    "r5", m2 → m4 + m6, 0, 100
    "r6", m4 → nothing, 0, 100
    "biomass", m6 + m5 → nothing, 0, 100
end

gs = [Gene("g$i") for i = 1:4]

model.reactions["biomass"].objective_coefficient = 1.0

add_genes!(model, gs)
add_metabolites!(model, [m1, m2, m3, m4, m5, m6])

reaction_isozymes = Dict(
    "r3" => [Isozyme(Dict("g1" => 1), 10.0, 10.0)],
    "r4" => [Isozyme(Dict("g2" => 1, "g3" => 3), 30.0, 20.0)],
    "r5" => [Isozyme(Dict("g3" => 1, "g4" => 2), 70.0, 30.0)],
)
gene_product_bounds = Dict(
    "g1" => (0.0, 0.2),
    "g2" => (0.0, 0.1),
    "g3" => (0.0, 10.0),
    "g4" => (0.0, 1000.0),
)

gene_product_molar_mass = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0)

gene_product_mass_group_bound = Dict("uncategorized" => 1.0)

model

# Construct and simulate a GECKO model
gecko_model = make_gecko_model(
    model;
    reaction_isozymes,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound,
)

# Get classic GECKO solution
optimized_model = flux_balance_analysis(
    gecko_model,
    CPLEX.Optimizer;
)
gecko_fluxes = flux_dict(gecko_model, optimized_model)
gecko_gps = gene_product_dict(gecko_model, optimized_model)
gene_product_mass_group_dict(gecko_model, optimized_model)

# Differentiate an optimal solution
rid_enzyme = Dict(
    k => isozyme_to_enzyme(first(v), gene_product_molar_mass; direction = :forward)
    for (k, v) in reaction_isozymes
)

diffmodel = with_parameters(gecko_model, rid_enzyme)

x, dx = differentiate(
    diffmodel,
    CPLEX.Optimizer
)
```

# More complicated examples
