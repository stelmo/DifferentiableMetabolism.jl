
# # Differentiating enzyme constrained metabolic models 

using DifferentiableMetabolism

using Symbolics
using ConstraintTrees
using COBREXA
using Tulip

include(joinpath("..", "test", "simple_model.jl")) #hide

# ![simple_model](./assets/simple_model.svg)

model

# ## Add enzyme kinetic information

Symbolics.@variables kcats_forward[1:4] kcats_backward[1:4]

reaction_isozymes = Dict(
    "r3" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[1], kcats_backward[1]),
    ),
    "r4" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[2], kcats_backward[2]),
        "iso2" =>
            ParameterIsozyme(Dict("g2" => 1), kcats_forward[3], kcats_backward[3]),
    ),
    "r5" => Dict(
        "iso1" => ParameterIsozyme(
            Dict("g3" => 1, "g4" => 2),
            kcats_forward[4],
            kcats_backward[4],
        ),
    ),
)

gene_molar_masses = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0, "g5" => 1.0)

 # ## Add a parameterized capacity limitation
 
Symbolics.@variables capacitylimitation

