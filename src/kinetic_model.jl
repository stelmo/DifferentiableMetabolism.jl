function build_kinetic_model(
    model::AbstractFBCModels.AbstractFBCModel;
    reaction_isozymes::Union{
        Dict{String,Dict{String,Isozyme}},
        Dict{String,Dict{String,ParameterIsozyme}},
    },
    gene_product_molar_masses::Dict{String,Float64},
    capacity::Union{Vector{Tuple{String,Vector{String},Real}},Real},
)
    # prepare some accessor functions for the later stuff
    # TODO: might be nicer to somehow parametrize the fwd/rev directions out.
    # Also there is a lot of conversion between symbols and strings, might be
    # nicer to have that sorted out in some better way.
    function isozyme_forward_ids(rid)
        haskey(reaction_isozymes, String(rid)) || return nothing
        return [
            Symbol(k) for
            (k, i) in reaction_isozymes[String(rid)] if !isnothing(i.kcat_forward)
        ]
    end
    function isozyme_reverse_ids(rid)
        haskey(reaction_isozymes, String(rid)) || return nothing
        return [
            Symbol(k) for
            (k, i) in reaction_isozymes[String(rid)] if !isnothing(i.kcat_reverse)
        ]
    end
    kcat_forward(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_forward
    kcat_reverse(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_reverse
    isozyme_gene_product_stoichiometry(rid, iso_id) = Dict(
        Symbol(k) => v for (k, v) in
        reaction_isozymes[String(rid)][String(iso_id)].gene_product_stoichiometry
    )
    gene_ids = Symbol.(keys(gene_product_molar_masses))
    gene_product_molar_mass(gid) = get(gene_product_molar_masses, String(gid), 0.0)

    # allocate all variables and build the system
    constraints = flux_balance_constraints(model)

    constraints += sign_split_variables(
        constraints.fluxes,
        positive = :fluxes_forward,
        negative = :fluxes_reverse,
    )

    constraints += enzyme_variables(;
        fluxes_forward = constraints.fluxes_forward,
        fluxes_reverse = constraints.fluxes_reverse,
        isozyme_forward_ids,
        isozyme_reverse_ids,
    )

    return constraints *
           :directional_flux_balance^sign_split_constraints(;
               positive = constraints.fluxes_forward,
               negative = constraints.fluxes_reverse,
               signed = constraints.fluxes,
           ) *
           enzyme_constraints(;
               fluxes_forward = constraints.fluxes_forward,
               fluxes_reverse = constraints.fluxes_reverse,
               isozyme_forward_amounts = constraints.isozyme_forward_amounts,
               isozyme_reverse_amounts = constraints.isozyme_reverse_amounts,
               kcat_forward,
               kcat_reverse,
               isozyme_gene_product_stoichiometry,
               gene_product_molar_mass,
               capacity_limits = capacity isa Real ?
                                 [(
                   :total_capacity,
                   gene_ids,
                   capacity isa Float64 ? C.Between(0, capacity) :
                   ParameterBetween(0, capacity),
               )] :
                                 [
                   (
                       Symbol(k),
                       Symbol.(gs),
                       capacity isa Float64 ? C.Between(0, capacity) :
                       ParameterBetween(0, cap),
                   ) for (k, gs, cap) in capacity
               ],
           )
end

export build_kinetic_model
