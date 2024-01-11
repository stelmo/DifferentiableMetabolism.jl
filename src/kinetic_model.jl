
function build_kinetic_model(
    model::AbstractFBCModels.AbstractFBCModel;
    reaction_isozymes::Union{
        Dict{String,Dict{String,COBREXA.Isozyme}},
        Dict{String,Dict{String,ParameterIsozyme}},
    },
    gene_product_molar_masses::Union{Dict{String,Float64},Dict{String,Symbolics.Num}},
    capacity::Union{
        Vector{Tuple{String,Vector{String},Float64}},
        Float64,
        Symbolics.Num,
        Vector{Tuple{String,Vector{String},Symbolics.Num}},
    },
)
    constraints = COBREXA.fbc_model_constraints(model)

    # might be nice to omit some conditionally (e.g. slash the direction if one
    # kcat is nothing)
    isozyme_amounts = COBREXA.isozyme_amount_variables(
        Symbol.(keys(reaction_isozymes)),
        rid -> Symbol.(keys(reaction_isozymes[string(rid)])),
    )

    # allocate variables for everything (nb. += wouldn't associate right here)
    constraints =
        constraints +
        :fluxes_forward^COBREXA.unsigned_positive_contribution_variables(
            constraints.fluxes,
        ) +
        :fluxes_reverse^COBREXA.unsigned_negative_contribution_variables(
            constraints.fluxes,
        ) +
        :isozyme_forward_amounts^isozyme_amounts +
        :isozyme_reverse_amounts^isozyme_amounts +
        :gene_product_amounts^ConstraintTrees.variables(
            keys = Symbol.(AbstractFBCModels.genes(model)),
            bounds = ConstraintTrees.Between(0, Inf),
        )

    # connect all parts with constraints
    constraints =
        constraints *
        :directional_flux_balance^COBREXA.sign_split_constraints(
            positive = constraints.fluxes_forward,
            negative = constraints.fluxes_reverse,
            signed = constraints.fluxes,
        ) *
        :isozyme_flux_forward_balance^COBREXA.isozyme_flux_constraints(
            constraints.isozyme_forward_amounts,
            constraints.fluxes_forward,
            (rid, isozyme) -> COBREXA.maybemap(
                x -> x.kcat_forward,
                COBREXA.maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
        ) *
        :isozyme_flux_reverse_balance^COBREXA.isozyme_flux_constraints(
            constraints.isozyme_reverse_amounts,
            constraints.fluxes_reverse,
            (rid, isozyme) -> COBREXA.maybemap(
                x -> x.kcat_reverse,
                COBREXA.maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
        ) *
        :gene_product_isozyme_balance^COBREXA.gene_product_isozyme_constraints(
            constraints.gene_product_amounts,
            (constraints.isozyme_forward_amounts, constraints.isozyme_reverse_amounts),
            (rid, isozyme) -> COBREXA.maybemap(
                x -> [(Symbol(k), v) for (k, v) in x.gene_product_stoichiometry],
                COBREXA.maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
        ) *
        :gene_product_capacity^(
            capacity isa Real ?
            ConstraintTrees.Constraint(
                value = sum(
                    gpa.value * gene_product_molar_masses[String(gp)] for
                    (gp, gpa) in constraints.gene_product_amounts
                ),
                bound = capacity isa Symbolics.Num ? ParameterBetween(0, capacity) :
                        ConstraintTrees.Between(0, capacity),
            ) :
            C.ConstraintTree(
                Symbol(id) => ConstraintTrees.Constraint(
                    value = sum(
                        constraints.gene_product_amounts[Symbol(gp)].value *
                        gene_product_molar_masses[gp] for gp in gps
                    ),
                    bound = limit isa Symbolics.Num ? ParameterBetween(0, limit) :
                            ConstraintTrees.Between(0, limit),
                ) for (id, gps, limit) in capacity_limits
            )
        )

    constraints
end

export build_kinetic_model
