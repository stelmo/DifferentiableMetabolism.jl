# DifferentiableMetabolism.jl
This package extends [COBREXA.jl](https://github.com/LCSB-BioCore/COBREXA.jl)
with the ability to differentiate an optimal solution of a constraint-based
metabolic model with respect to parameters. 

However, please note that only non-degenerate solutions (unique) can be
differentiated for the derivatives to have a concrete interpretation. For enzyme
constrained models, this means that you will only be able to differentiate the
model if you prune the inactive reactions from the solution (see the
documentation for examples). For other types of models, you need to ensure that
your solution is non-degenerate, otherwise you will only compute sub-gradients.

To use this package you

To use this package, [download and install Julia](https://julialang.org/downloads/), and add 
the following packages using the built in package manager:
```julia
] add COBREXA, DifferentiableMetabolism, Tulip
```
Note, any optimization solver that is compatible with [JuMP](https://jump.dev/)
can be used. Here we have opted to use
[Tulip.jl](https://github.com/ds4dm/Tulip.jl). To run the tests,
[Ipopt](https://github.com/jump-dev/Ipopt.jl) is also required.
```julia
] test DifferentiableMetabolism
```
## Differentiating a simple model
In this example, a simple model will be differentiated.

```julia
using DifferentiableMetabolism
using Tulip
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
    "r7", m2 → m4 + m6, 0, 100
    "biomass", m6 + m5 → nothing, 0, 100
end

gs = [Gene("g$i") for i = 1:5]

model.reactions["biomass"].objective_coefficient = 1.0

add_genes!(model, gs)
add_metabolites!(model, [m1, m2, m3, m4, m5, m6])

reaction_isozymes = Dict(
    "r3" => [Isozyme(Dict("g1" => 1), 10.0, 10.0)],
    "r4" => [Isozyme(Dict("g2" => 1, "g3" => 3), 30.0, 20.0)],
    "r5" => [Isozyme(Dict("g3" => 1, "g4" => 2), 70.0, 30.0)],
    "r7" => [Isozyme(Dict("g5" => 1), 50.0, 20.0)],
)
gene_product_bounds = Dict(
    "g1" => (0.0, 0.2),
    "g2" => (0.0, 0.1),
    "g3" => (0.0, 10.0),
    "g4" => (0.0, 1000.0),
    "g5" => (0.0, 1000.0),
)

gene_product_molar_mass = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0, "g5" => 5.0)

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
    Tulip.Optimizer;
)
gecko_fluxes = flux_dict(gecko_model, optimized_model) # notice that r5 is inactive!

# Prune away inactive reactions
pruned_model = prune_model(model, gecko_fluxes)

# Differentiate an optimal solution
pruned_gecko_model = make_gecko_model(
    pruned_model;
    reaction_isozymes,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound,
)

rid_enzyme = Dict(
    k => isozyme_to_enzyme(first(v), gene_product_molar_mass; direction = :forward)
    for (k, v) in reaction_isozymes
)

diffmodel = with_parameters(gecko_model, rid_enzyme)

x, dx = differentiate(
    diffmodel,
    Tulip.Optimizer
)
```
Here, `x` are the variables, corresponding to `diffmodel.var_ids`, and `dx` are
the derivatives, where rows correspond to `diffmodel.param_ids`, and columns
correspond to `diffmodel.var_ids`.

While this package is under development, you can already use more advanced
functionality. Look at the tests to see how to incorporate thermodynamic and/or 
saturation effects these differentiable models.
