"""
    _remove_lin_dep_rows(A; ϵ = 1e-8, atol = 1e-8, digits = 16)

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
    _dg(
        model,
        rid_enzyme,
        rid_dg0,
        rid,
        θ;
        RT = 298.15 * 8.314e-3,
        ignore_reaction_ids = [],
    )

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
    _saturation(
        model,
        rid_enzyme,
        rid_km,
        rid,
        mangled_rid,
        θ;
        ignore_reaction_ids = [],
    )

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
    prune_model(
        model::StandardModel,
        reaction_fluxes,
        gene_product_concentrations = Dict{String,Float64}();
        atol = 1e-9,
        verbose = true,
    )

Return a simplified version of `model` that contains only reactions (and the
associated metabolites) that are active, i.e. carry fluxes (from
`reaction_fluxes`) absolutely bigger than `atol`. All reactions are set
unidirectional based on `reaction_fluxes`. If `gene_product_concentrations` is
supplied, then genes that have a concentration bigger than `atol` are also
included in the pruned model, otherwise no gene information is retained.
"""
function prune_model(
    model::StandardModel,
    reaction_fluxes,
    gene_product_concentrations = Dict{String,Float64}();
    atol = 1e-9,
    verbose = true,
)
    pruned_model = StandardModel("pruned_model")

    rxns = Vector{Reaction}()
    mets = Vector{Metabolite}()
    gs = Vector{Gene}()
    mids = String[]

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
    end

    for mid in unique(mids)
        push!(mets, model.metabolites[mid])
    end

    for (gid, conc) in gene_product_concentrations
        abs(conc) <= atol && continue
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
    _make_differentiable_model(
        c,
        _E,
        d,
        M,
        h,
        θ,
        var_ids,
        param_ids;
        scale_equality = false,
    )

Internal helper function to construct a basic linear program representing a
differentiable metabolic model.
"""
function _make_differentiable_model(
    c,
    _E,
    d,
    M,
    h,
    θ,
    var_ids,
    param_ids;
    scale_equality = false,
)

    Q = _ -> spzeros(length(var_ids), length(var_ids))

    if scale_equality
        row_factors = scaling_factor(_E(θ), d(θ))
    else
        row_factors = fill(1.0, size(_E(θ), 1))
    end

    E(θ) = row_factors .* _E(θ)

    return DifferentiableModel(
        Q,
        c,
        E,
        d,
        M,
        h,
        θ,
        _ -> throw(MissingException("Missing method: analytic derivatives of variables.")),
        _ -> throw(MissingException("Missing method: analytic derivatives of parameters.")),
        var_ids,
        param_ids,
    )
end
