"""
"""
function differentiate_thermodynamic_smoment(
    smm::SMomentModel,
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

    kcat_rid_order = [
        rid for
        rid in reactions(smm.smodel) if haskey(smm.enzymedata.reaction_kcats, rid) &&
        COBREXA.has_reaction_grr(smm.smodel, rid)
    ]

    met_conc_order =
        [mid for mid in metabolites(smm.smodel) if haskey(metabolite_concentrations, mid)]

    θ = [
        [first(smm.enzymedata.reaction_kcats[rid][1]) for rid in kcat_rid_order]
        [metabolite_concentrations[mid] for mid in met_conc_order]
        smm.enzymedata.total_protein_mass
    ]

    opt, opt_struct = differentiable_thermodynamic_smoment_opt_problem(
        smm;
        reaction_dg0s,
        kcat_rid_order,
        met_conc_order,
        ϵ,
        RT,
        ignore_reaction_ids,
        stoich_digits_round,
    )

    if isnothing(objective_id)
        obj_idx_orig = first(findnz(objective(smm.smodel))[1])
        obj_id_orig = reactions(smm.smodel)[obj_idx_orig]
        obj_idx = opt_struct.reaction_map[obj_id_orig*"§FOR"]
    else
        obj_idx = opt_struct.reaction_map[objective_id*"§FOR"]
    end
    opt.c[obj_idx] = -1.0 # max due to min sense

    x, dx, _ = differentiate_LP(
        opt,
        θ,
        optimizer;
        use_analytic = false, # no analytic option
        scale_input,
        scale_output,
        modifications,
    )

    # update smoment internals
    smm.smomentdata = SMomentData(
        opt.c,
        opt.Ef(θ),
        opt.d,
        opt.M,
        opt.hf(θ),
        opt_struct.reaction_map,
        opt_struct.metabolite_map,
    )

    return (
        x = x,
        dx = dx,
        θ = θ,
        rids = COBREXA._order_id_to_idx_dict(opt_struct.reaction_map),
        kcat_rid_order,
        met_conc_order,
    )
end

"""
Return structures that will allow the most basic form of 
smoment to be solved. No enzyme constraints allowed.
Most effective enzyme is the only GRR. Assume unidirectional reactions.
"""
function differentiable_thermodynamic_smoment_opt_problem(
    smm::SMomentModel;
    reaction_dg0s = Dict(),
    kcat_rid_order = [],
    met_conc_order = [],
    ϵ = 1e-8,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
    stoich_digits_round = 8,
)

    S_rank_issues, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        COBREXA._build_irreversible_stoichiometric_matrix(smm.smodel)

    #: make full rank
    S = round.(_remove_lin_dep_rows(S_rank_issues; ϵ), digits = stoich_digits_round)

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
        !haskey(smm.enzymedata.reaction_kcats, original_rid) && continue
        original_rid in ignore_reaction_ids && continue

        # these entries have kcats, only one GRR by assumption
        grr = first(reaction_gene_association(smm.smodel, original_rid))
        pstoich = first(smm.enzymedata.reaction_protein_stoichiometry[original_rid])
        mw = dot(pstoich, [smm.enzymedata.protein_masses[gid] for gid in grr])
        rs = reaction_stoichiometry(smm.smodel, original_rid)
        dir = COBREXA.is_reaction_forward_only(smm.smodel, original_rid) ? 1.0 : -1.0 # unidirectional by assumption

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
    d = spzeros(n_metabolites + 1)

    #: need to set objective reaction outside
    c = spzeros(n_vars)

    #: inequality constraints
    M, h = COBREXA._smoment_build_inequality_constraints(
        smm,
        n_reactions,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )
    hf(θ) = [h[1:end-1]; last(θ)]

    opt = (c = c, Ef = Ef, d = d, M = sparse(M), hf = hf)
    opt_struct = (
        reaction_map = reaction_map,
        metabolite_map = metabolite_map,
        n_reactions = n_reactions,
        n_metabolites = n_metabolites,
        kcat_rid_order = kcat_rid_order,
        met_conc_order = met_conc_order,
    )

    return opt, opt_struct
end
