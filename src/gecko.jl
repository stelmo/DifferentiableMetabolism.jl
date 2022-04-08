"""
    differentiate_gecko(
        gm::GeckoModel,
        optimizer;
        objective_id = nothing,
        ϵ = 1e-8,
        stoich_digits_round = 8,
        use_analytic = false,
        scale_input = false,
        scale_output = true,
        modifications = [],
    )

Notes:
1) `model` must be pruned, i.e. only active reactions, metabolites and genes.
Singularity issues will arise if not.
2) Each reaction has only one active isozyme and one kcat assigned to that
isozyme.
3) Rescale some columns to prevent numerical issues from arising
"""
function differentiate_gecko(
    gm::GeckoModel,
    optimizer;
    objective_id = nothing,
    ϵ = 1e-8,
    stoich_digits_round = 8,
    use_analytic = false,
    scale_input = false,
    scale_output = true,
    modifications = [],
)

    kcat_rid_order = [
        rid for
        rid in reactions(gm.smodel) if haskey(gm.enzymedata.reaction_kcats, rid) &&
        COBREXA.has_reaction_grr(gm.smodel, rid)
    ]
    θ = [
        [first(gm.enzymedata.reaction_kcats[rid][1]) for rid in kcat_rid_order]
        gm.enzymedata.total_protein_mass
    ]

    opt, opt_struct =
        differentiable_gecko_opt_problem(gm; kcat_rid_order, ϵ, stoich_digits_round)

    if isnothing(objective_id)
        obj_idx_orig = first(findnz(objective(gm.smodel))[1])
        obj_id_orig = reactions(gm.smodel)[obj_idx_orig]
        obj_idx = opt_struct.reaction_map[obj_id_orig*"§FOR"]
    else
        obj_idx = opt_struct.reaction_map[objective_id*"§FOR"]
    end
    opt.c[obj_idx] = -1.0 # max due to min sense

    if use_analytic
        dFdθ = analytic_gecko_derivatives(opt_struct)
    else
        dFdθ = nothing
    end
    x, dx, _ = differentiate_LP(
        opt,
        θ,
        optimizer;
        use_analytic,
        scale_input,
        scale_output,
        modifications,
        dFdθ,
    )

    # update geckomodel internals
    gm.geckodata = GeckoData(
        opt.c,
        opt.Ef(θ),
        opt.d,
        opt.M,
        opt.hf(θ),
        opt_struct.reaction_map,
        opt_struct.metabolite_map,
        opt_struct.protein_ids,
    )

    return (
        x = x,
        dx = dx,
        θ = θ,
        rids = COBREXA._order_id_to_idx_dict(gm.geckodata.reaction_map),
        pids = gm.geckodata.protein_ids,
    )
end

"""
Return optimization problem for gecko problem, but in differentiable format.
"""
function differentiable_gecko_opt_problem(
    model::GeckoModel;
    kcat_rid_order = [],
    ϵ = 1e-8,
    stoich_digits_round = 8,
)
    S_rank_issues, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        COBREXA._build_irreversible_stoichiometric_matrix(model.smodel)

    #: make full rank
    S = round.(_remove_lin_dep_rows(S_rank_issues; ϵ), digits = stoich_digits_round)

    #: find all gene products that have kcats associated with them
    protein_ids = COBREXA._get_proteins_with_kcats(model)

    #: size of resultant model
    n_reactions = size(S, 2)
    n_proteins = length(protein_ids)
    n_metabolites = size(S, 1)
    n_vars = n_reactions + n_proteins

    #: equality lhs
    E_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeff_tuple = Vector{Tuple{Float64,String}}(),
    )

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        if contains(rid, "§ARM")
            @warn "Isozyme found, unexpected!"
            continue
        end
        !haskey(model.enzymedata.reaction_kcats, original_rid) && continue

        # add all entries to column of matrix
        #! no isozymes by assumption and all unidirectional
        _add_enzyme_variable_as_function(
            model,
            original_rid,
            E_components,
            col_idx,
            protein_ids,
        )
    end

    kcat_pst_idxs = [
        (pst, first(indexin([rid], kcat_rid_order))) for
        (pst, rid) in E_components.coeff_tuple
    ]

    Se(θ) = sparse(
        E_components.row_idxs,
        E_components.col_idxs,
        [pst / θ[x] for (pst, x) in kcat_pst_idxs],
        n_proteins,
        n_reactions,
    )

    Ef(θ) = [
        S spzeros(n_metabolites, n_proteins)
        Se(θ) I(n_proteins)
    ]

    #: equality rhs
    d = spzeros(n_metabolites + n_proteins) #TODO handle fixed variables

    #: need to set objective reaction outside
    c = spzeros(n_vars)

    #: inequality constraints
    M, h = COBREXA._gecko_build_inequality_constraints(
        model,
        protein_ids,
        n_reactions,
        n_proteins,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )
    hf(θ) = [h[1:end-1]; last(θ)]

    opt = (c = c, Ef = Ef, d = d, M = sparse(M), hf = hf)
    opt_struct = (
        reaction_map = reaction_map,
        metabolite_map = metabolite_map,
        protein_ids = protein_ids,
        kcat_rid_order = kcat_rid_order,
        row_idxs = E_components.row_idxs,
        col_idxs = E_components.col_idxs,
        vals = kcat_pst_idxs,
        n_proteins = n_proteins,
        n_reactions = n_reactions,
        n_metabolites = n_metabolites,
    )

    return opt, opt_struct
end

"""
Helper function to add an column into the enzyme stoichiometric matrix
parametrically. Note, assumes only one grr.
"""
function _add_enzyme_variable_as_function(
    model::GeckoModel,
    original_rid,
    E_components,
    col_idx,
    protein_ids,
)
    grr = reaction_gene_association(model, original_rid)[1] #! assumes both kcats are the same (appropriate one duplicated)
    pstoich = model.enzymedata.reaction_protein_stoichiometry[original_rid][1] #! requires pruned protein_stoichiometry
    for (idx, pid) in enumerate(grr)
        push!(E_components.row_idxs, first(indexin([pid], protein_ids)))
        push!(E_components.col_idxs, col_idx)
        push!(E_components.coeff_tuple, (-pstoich[idx], original_rid))
    end
end

"""
Helper function to get gecko problem derivatives analytically.
"""
function _block1_gecko_derivatives(z, ν, λ, θ, opt_struct)
    row_idxs = Int[]
    col_idxs = Int[]
    vals = Float64[]
    ν_star = ν[(1+opt_struct.n_metabolites):end]

    for (idx, jdx, v) in zip(opt_struct.row_idxs, opt_struct.col_idxs, opt_struct.vals)
        coeff, kcat_idx = v
        push!(col_idxs, kcat_idx)
        push!(row_idxs, jdx)
        push!(vals, ν_star[idx] * coeff * 1.0 / θ[kcat_idx]^2) # negative multiplied through
    end
    n_rows = opt_struct.n_reactions
    n_cols = length(θ)
    sparse(row_idxs, col_idxs, vals, n_rows, n_cols)
end

"""
Helper function to get gecko problem derivatives analytically.
"""
function _block2_gecko_derivatives(z, ν, λ, θ, opt_struct)
    row_idxs = Int[]
    col_idxs = Int[]
    vals = Float64[]

    for (idx, jdx, v) in zip(opt_struct.row_idxs, opt_struct.col_idxs, opt_struct.vals)
        coeff, kcat_idx = v
        push!(col_idxs, kcat_idx)
        push!(row_idxs, idx)
        push!(vals, z[jdx] * coeff * -1.0 / θ[kcat_idx]^2)
    end
    n_rows = opt_struct.n_proteins
    n_cols = length(θ)
    sparse(row_idxs, col_idxs, vals, n_rows, n_cols)
end

"""
Analytic derivative of KKT conditions for GECKO type problems. 
Much faster than autodiff approach.
"""
function analytic_gecko_derivatives(opt_struct)
    (z, ν, λ, θ) -> Array(
        [
            _block1_gecko_derivatives(z, ν, λ, θ, opt_struct)
            zeros(opt_struct.n_proteins, length(θ))
            zeros(opt_struct.n_metabolites, length(θ))
            _block2_gecko_derivatives(z, ν, λ, θ, opt_struct)
            zeros(length(λ) - 1, length(θ))
            [zeros(length(θ) - 1); -λ[end]]'
        ],
    )
end
