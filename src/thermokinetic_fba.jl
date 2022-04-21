function differentiate_thermokinetic(
    model::StandardModel,
    optimizer;
    reaction_dg0s = Dict(),
    metabolite_concentrations = Dict(),
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
    modifications = [],
    ϵ = 1e-8,
    objective_id = nothing,
    stoich_digits_round = 8,
    scale_input = false,
    scale_output = true,
)

end
