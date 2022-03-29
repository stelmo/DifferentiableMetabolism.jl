function gecko_opt_problem_temp(
    model::StandardModel;
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    lb_protein_measurements = Dict(),
    ub_protein_measurements = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    total_protein_mass = 0.0,
)
    S, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        COBREXA._build_irreversible_stoichiometric_matrix(model)

    #: find all gene products that have kcats associated with them
    protein_ids = COBREXA._get_proteins_with_kcats(model, reaction_kcats)

    #: size of resultant model
    n_reactions = size(S, 2)
    n_proteins = length(protein_ids)
    n_metabolites = size(S, 1)
    n_vars = n_reactions + n_proteins

    #: equality lhs
    E_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeffs = Vector{Float64}(),
    )

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        contains(rid, "§ARM") && continue
        !haskey(reaction_kcats, original_rid) && continue

        # these entries have kcats
        if contains(rid, "§ISO")
            iso_num = parse(
                Int,
                replace(
                    first(filter(startswith("ISO"), split(rid, "§")[2:end])),
                    "ISO" => "",
                ),
            )
        else # only one enzyme
            iso_num = 1
        end

        # add all entries to column of matrix
        COBREXA._add_enzyme_variable(
            model,
            iso_num, # only one enzyme
            rid,
            original_rid,
            protein_stoichiometry,
            reaction_kcats,
            E_components,
            col_idx,
            protein_ids,
        )
    end

    Se = sparse(
        E_components.row_idxs,
        E_components.col_idxs,
        E_components.coeffs,
        n_proteins,
        n_reactions,
    )

    E = [
        S zeros(n_metabolites, n_proteins)
        Se I(n_proteins)
    ]

    #: equality rhs
    d = zeros(n_metabolites + n_proteins)

    #: need to set objective reaction outside
    c = spzeros(n_vars)

    #: inequality constraints
    M, h = COBREXA._gecko_build_inequality_constraints(
        lb_protein_measurements,
        ub_protein_measurements,
        protein_ids,
        protein_masses,
        n_reactions,
        n_proteins,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
        total_protein_mass,
    )

    return c, E, d, M, h, reaction_map, metabolite_map, protein_ids
end
