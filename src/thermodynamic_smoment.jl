"""
"""
function thermodynamic_smoment(
    model::StandardModel,
    optimizer;
    objective_id = "",
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    reaction_dg0s = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    metabolite_concentrations = Dict(),
    total_protein_mass = 0.0,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
    sense = MOI.MAX_SENSE,
    modifications = [],
)
    _, E, d, M, h, reaction_map, _ = thermodynamic_smoment_opt_problem(
        model;
        protein_stoichiometry,
        protein_masses,
        reaction_kcats,
        reaction_dg0s,
        lb_flux_measurements,
        ub_flux_measurements,
        metabolite_concentrations,
        total_protein_mass,
        RT,
        ignore_reaction_ids,
    )

    opt_model = Model(optimizer)
    x = @variable(opt_model, x[1:size(E, 2)])
    bid = reaction_map[objective_id]
    @objective(opt_model, sense, x[bid])
    @constraint(opt_model, E * x .== d)
    @constraint(opt_model, M * x .<= h)

    # apply the modifications, if any
    for mod in modifications
        mod(nothing, opt_model)
    end

    optimize!(opt_model)

    COBREXA._map_irrev_to_rev_ids(reaction_map, value.(x))
end

"""

Return structures that will allow the most basic form of 
thermokinetic smoment to be solved. No enzyme constraints allowed.
Most effective enzyme is the only GRR.

Assume unidirectional reactions!
"""
function thermodynamic_smoment_opt_problem(
    model::StandardModel;
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    reaction_dg0s = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    metabolite_concentrations = Dict(),
    total_protein_mass = 0.0,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
)

    S, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        COBREXA._build_irreversible_stoichiometric_matrix(model)

    #: size of resultant model
    n_reactions = size(S, 2)
    n_metabolites = size(S, 1)
    n_vars = n_reactions + 1

    #: equality lhs
    Se = zeros(1, n_reactions)

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        !haskey(reaction_kcats, original_rid) && continue
        original_rid in ignore_reaction_ids && continue

        # these entries have kcats, only one GRR by assumption
        grr = first(reaction_gene_association(model, original_rid))
        pstoich = first(protein_stoichiometry[original_rid])
        mw = dot(pstoich, [protein_masses[gid] for gid in grr])
        kcat =
            contains(rid, "§FOR") ? first(reaction_kcats[original_rid])[1] :
            first(reaction_kcats[original_rid])[2]

        thermo_factor = 1.0
        if haskey(reaction_dg0s, original_rid)
            rstoich = reaction_stoichiometry(model, original_rid)
            if all(in.(keys(rstoich), Ref(keys(metabolite_concentrations))))
                dg =
                    reaction_dg0s[original_rid] +
                    RT * sum(
                        rstoich[mid] * log(metabolite_concentrations[mid]) for
                        mid in keys(rstoich)
                    )
                sf = dg < 0.0 ? 1.0 : -1.0
                thermo_factor = 1 - exp(sf * dg / RT)
            end
        end
        Se[1, col_idx] = -mw / (kcat * thermo_factor)
    end

    E = [
        S zeros(n_metabolites, 1)
        Se 1.0
    ]

    #: equality rhs
    d = zeros(n_metabolites + 1)

    #: need to set objective reaction outside
    c = spzeros(n_vars)

    #: inequality constraints
    M, hf = _smoment_build_inequality_constraints_as_functions(
        n_reactions,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

    return c, E, d, M, hf([total_protein_mass]), reaction_map, metabolite_map
end

"""
Return structures that will allow the most basic form of 
smoment to be solved. No enzyme constraints allowed.
Most effective enzyme is the only GRR. Assume unidirectional reactions.
"""
function differentiable_thermokinetic_smoment(
    model::StandardModel;
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_dg0s = Dict(),
    reaction_kcats = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    kcat_rid_order = [],
    met_conc_order = [],
    ϵ = 1e-8,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
)

    S_rank_issues, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        COBREXA._build_irreversible_stoichiometric_matrix(model)

    #: make full rank
    S = _remove_lin_dep_rows(S_rank_issues; ϵ)

    #: size of resultant model
    n_reactions = size(S, 2)
    n_metabolites = size(S, 1)
    n_vars = n_reactions + 1

    #: equality lhs
    kcat_original_rid_order = String[]
    met_order = Vector{Vector{String}}()
    nus_vec = Vector{Vector{Float64}}()
    col_idxs = Int[]
    neg_mws = Float64[]
    dirs = Float64[]
    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        !haskey(reaction_kcats, original_rid) && continue
        original_rid in ignore_reaction_ids && continue

        # these entries have kcats, only one GRR by assumption
        grr = first(reaction_gene_association(model, original_rid))
        pstoich = first(protein_stoichiometry[original_rid])
        mw = dot(pstoich, [protein_masses[gid] for gid in grr])
        rs = reaction_stoichiometry(model, original_rid)
        dir = COBREXA._is_forward_only(model, original_rid) ? 1.0 : -1.0 # unidirectional by assumption

        # kinetic info
        push!(kcat_original_rid_order, original_rid)
        push!(col_idxs, col_idx)
        push!(neg_mws, -mw)

        # thermo info
        if haskey(reaction_dg0s, original_rid)
            push!(met_order, collect(keys(rs)))
            push!(nus_vec, collect(values(rs)))
            push!(dirs, dir)
        else
            push!(met_order, collect(keys(rs)))
            push!(nus_vec, Float64[]) # acts like a flag
            push!(dirs, dir)
        end
    end

    kcat_idxs = [
        (first(indexin([rid], kcat_rid_order)), nmw) for
        (rid, nmw) in zip(kcat_original_rid_order, neg_mws)
    ]

    mid_idxs = [
        (dir, nus, Int.(indexin(mids, met_conc_order)) .+ length(kcat_rid_order)) for
        (dir, nus, mids) in zip(dirs, nus_vec, met_order)
    ]

    function dg(dir, rid_idx, nus, mconcs)
        isempty(nus) && return 1.0

        dg_val =
            reaction_dg0s[kcat_rid_order[rid_idx]] +
            RT * sum(nu * log(mconc) for (nu, mconc) in zip(nus, mconcs))
        return 1 - exp(dir * dg_val / RT)
    end

    fSe(θ) = sparsevec(
        col_idxs,
        [
            nmw / ((θ[rid_idx]) * dg(dir, rid_idx, nus, θ[mconcs_idxs])) for
            ((rid_idx, nmw), (dir, nus, mconcs_idxs)) in zip(kcat_idxs, mid_idxs)
        ],
        n_reactions,
    )

    Ef(θ) = [
        Array(S) zeros(n_metabolites, 1)
        fSe(θ)' 1.0
    ]

    #: equality rhs
    d = zeros(n_metabolites + 1)

    #: need to set objective reaction outside
    c = zeros(n_vars)

    #: inequality constraints
    M, hf = _smoment_build_inequality_constraints_as_functions(
        n_reactions,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

    return c, Ef, d, M, hf, reaction_map, metabolite_map
end
