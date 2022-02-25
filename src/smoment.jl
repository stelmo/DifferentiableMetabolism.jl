"""
Return structures that will allow the most basic form of smoment to be solved.
No enzyme constraints allowed. Assume preprocessing changes model such that most
effective enzyme is the only GRR.
"""
function differentiable_smoment_opt_problem(
    model::StandardModel;
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    kcat_rid_order = [],
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    ϵ = 1e-8,
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
    col_idxs = Int[]
    neg_mws = Float64[]
    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        !haskey(reaction_kcats, original_rid) && continue

        # these entries have kcats, only one GRR by assumption
        grr = first(reaction_gene_association(model, original_rid))
        pstoich = first(protein_stoichiometry[original_rid])
        mw = dot(pstoich, [protein_masses[gid] for gid in grr])
        push!(kcat_original_rid_order, original_rid)
        push!(col_idxs, col_idx)
        push!(neg_mws, -mw)
    end

    kcat_idxs = [
        (first(indexin([rid], kcat_rid_order)), nmw) for
        (rid, nmw) in zip(kcat_original_rid_order, neg_mws)
    ]

    fSe(θ) = sparsevec(col_idxs, [nmw / θ[x] for (x, nmw) in kcat_idxs], n_reactions)

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

"""
    _smoment_build_inequality_constraints_as_functions(
        n_reactions,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

Helper function to return functions describing the inequality 
constraints for smoment.
"""
function _smoment_build_inequality_constraints_as_functions(
    n_reactions,
    lb_flux_measurements,
    ub_flux_measurements,
    lb_fluxes,
    ub_fluxes,
    reaction_map,
)
    #: inequality lhs
    M = Array(
        [
            -I(n_reactions) zeros(n_reactions, 1)
            I(n_reactions) zeros(n_reactions, 1)
            zeros(1, n_reactions) 1
        ],
    )

    #: inequality rhs
    for original_rid in keys(lb_flux_measurements) # only constrain if measurement available
        lb = lb_flux_measurements[original_rid]
        ub = ub_flux_measurements[original_rid]
        rids = [rid for rid in keys(reaction_map) if startswith(rid, original_rid)]

        if lb > 0 # forward only
            for rid in rids
                contains(rid, "§REV") && (ub_fluxes[reaction_map[rid]] = 0.0)
                contains(rid, "§FOR") &&
                    (ub_fluxes[reaction_map[rid]] = ub; lb_fluxes[reaction_map[rid]] = lb)
            end
        elseif ub < 0 # reverse only
            for rid in rids
                contains(rid, "§FOR") && (ub_fluxes[reaction_map[rid]] = 0.0)
                contains(rid, "§REV") &&
                    (ub_fluxes[reaction_map[rid]] = -lb; lb_fluxes[reaction_map[rid]] = -ub)
            end
        else # measurement does not rule our reversibility
            for rid in rids
                contains(rid, "§FOR") &&
                    (ub_fluxes[reaction_map[rid]] = ub; lb_fluxes[reaction_map[rid]] = 0)
                contains(rid, "§REV") &&
                    (ub_fluxes[reaction_map[rid]] = -lb; lb_fluxes[reaction_map[rid]] = 0)
            end
        end
    end

    hf(θ) = Array([-lb_fluxes; ub_fluxes; θ[end]])

    return M, hf
end
