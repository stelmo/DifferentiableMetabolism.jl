"""
$(TYPEDSIGNATURES)

Remove linearly dependent rows in `A` deduced by reducing the matrix to row
echelon form. Adjust sensitivity to numerical issues with `ϵ`. After row echelon
reduction, all elements are rounded to `digits` after the decimal. Rows of all zeros
(absolute value of the each element in row ≤ `atol`) are removed. Beware, this
operation is expensive for very large matrices.
"""
function _remove_lin_dep_rows(A; ϵ = 1e-8, atol = 1e-8, digits = 16)
    #TODO this method is suboptimal and can be improved for numerical stability
    #TODO improve RowEchelon, SVD does not work due to column reordering
    rA = round.(rref!(copy(Array(A)), ϵ); digits)

    idxs = Int[]
    for i = 1:size(rA, 1)
        if !all(abs.(rA[i, :]) .<= atol) # remove rows of all zero
            push!(idxs, i)
        end
    end

    return rA[idxs, :]
end

"""
$(TYPEDSIGNATURES)

Helper function to assign the thermodynamic driving force to a reaction.
"""
function _dg(
    model,
    rid_enzyme,
    rid_dg0,
    rid,
    θ;
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
)
    rs = reaction_stoichiometry(model, rid)
    stoich = values(rs)
    mids = collect(keys(rs))
    midxs = Int.(indexin(mids, metabolites(model))) .+ length(rid_enzyme)
    if !haskey(rid_dg0, rid) || rid in ignore_reaction_ids
        # no kinetic info or should be ignore thermodynamically
        1.0
    else
        dg_val =
            rid_dg0[rid] + RT * sum(nu * log(θ[midx]) for (nu, midx) in zip(stoich, midxs))
        1.0 - exp(dg_val / RT)
    end
end

"""
$(TYPEDSIGNATURES)
    
A helper function to incorporate saturation effects.
"""
function _saturation(
    model,
    rid_enzyme,
    rid_km,
    rid,
    mangled_rid,
    θ;
    ignore_reaction_ids = [],
)
    rs = reaction_stoichiometry(model, rid)
    stoich = values(rs)
    mids = collect(keys(rs))
    midxs = Int.(indexin(mids, metabolites(model))) .+ length(rid_enzyme)
    is_forward = contains(mangled_rid, "#forward") ? true : false
    @assert(contains(mangled_rid, rid)) #TODO sanity check, remove
    if !haskey(rid_km, rid) || rid in ignore_reaction_ids
        1.0
    else
        s_term = prod(
            (θ[midx] / rid_km[rid][mid])^nu for
            (nu, midx, mid) in zip(stoich, midxs, mids) if nu > 0
        )
        p_term = prod(
            (θ[midx] / rid_km[rid][mid])^nu for
            (nu, midx, mid) in zip(stoich, midxs, mids) if nu < 0
        )
        is_forward ? s_term / (1.0 + s_term + p_term) : p_term / (1.0 + s_term + p_term)
    end
end

"""
$(TYPEDSIGNATURES)

Return a simplified version of `model` that contains only reactions (and the
associated genes and metabolites) that are active, i.e. carry fluxes (from
`reaction_fluxes`) absolutely bigger than `atol`. All reactions are set
unidirectional based on `reaction_fluxes`.
"""
function prune_model(
    model::StandardModel,
    reaction_fluxes;
    atol = 1e-9,
    verbose = true,
)
    pruned_model = StandardModel("pruned_model")

    rxns = Vector{Reaction}()
    mets = Vector{Metabolite}()
    gs = Vector{Gene}()
    mids = String[]
    gids = String[]

    for rid in reactions(model)
        abs(reaction_fluxes[rid]) <= atol && continue

        rxn = deepcopy(model.reactions[rid])
        if reaction_fluxes[rid] > 0
            rxn.lb = max(0, rxn.lb)
        else
            rxn.ub = min(0, rxn.ub)
        end
        push!(rxns, rxn)

        rs = reaction_stoichiometry(model, rid)
        for mid in keys(rs)
            push!(mids, mid)
        end
        grrs = reaction_gene_association(model, rid)
        isnothing(grrs) && continue
        for grr in grrs
            append!(gids, grr)
        end
    end

    for mid in unique(mids)
        push!(mets, model.metabolites[mid])
    end

    for gid in unique(gids)
        push!(gs, model.genes[gid])
    end

    add_reactions!(pruned_model, rxns)
    add_metabolites!(pruned_model, mets)
    add_genes!(pruned_model, gs)

    #: print some info about process
    if verbose
        rrn = n_reactions(model) - n_reactions(pruned_model)
        @info "Removed $rrn reactions."
    end

    return pruned_model
end

"""
$(TYPEDSIGNATURES)

Internal helper function to construct a basic linear program representing a
differentiable metabolic model.
"""
function _make_differentiable_model(
    c,
    _E,
    _d,
    _M,
    _h,
    θ,
    var_ids,
    param_ids;
    scale_equality = false,
    scale_inequality = false,
)

    Q = _ -> spzeros(length(var_ids), length(var_ids))

    if scale_equality
        eq_row_factors = scaling_factor(_E(θ), _d(θ))
    else
        eq_row_factors = fill(1.0, size(_E(θ), 1))
    end

    if scale_inequality
        ineq_row_factors = scaling_factor(_M(θ), _h(θ))
    else
        ineq_row_factors = fill(1.0, size(_M(θ), 1))
    end

    E(θ) = eq_row_factors .* _E(θ)
    d(θ) = eq_row_factors .* _d(θ)

    M(θ) = ineq_row_factors .* _M(θ)
    h(θ) = ineq_row_factors .* _h(θ)

    return DifferentiableModel(
        Q,
        c,
        E,
        d,
        M,
        h,
        θ,
        (_, _, _, _) ->
            throw(MissingException("Missing method: analytic derivatives of variables.")),
        (_, _, _, _) ->
            throw(MissingException("Missing method: analytic derivatives of parameters.")),
        var_ids,
        param_ids,
    )
end
