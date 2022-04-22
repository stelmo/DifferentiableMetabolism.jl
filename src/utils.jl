"""
    $(TYPEDSIGNATURES)

Remove linearly dependent rows in `A` deduced by reducing the matrix to row
echelon form. Adjust sensitivity to numerical issues with `ϵ`. After row echelon
reduction, all elements are rounded to `digits` after the decimal. Rows of all zeros
(absolute value of the each element in row ≤ `atol`) are removed. Beware, this
operation is expensive for very large matrices.
"""
function _remove_lin_dep_rows(A; ϵ = 1e-8, atol = 1e-12, digits = 16)
    #TODO this method is suboptimal and can be improved for numerical stability
    #TODO improve RowEchelon, SVD does not work due to column reordering
    rA = round.(rref!(copy(Array(A)), ϵ); digits)

    idxs = Int[]
    for i = 1:size(rA, 1)
        if !all(abs.(rA[i, :]) .> atol) # remove rows of all zero
            push!(idxs, i)
        end
    end

    return rA[idxs, :]
end

"""
    $(TYPEDSIGNATURES)

Helper function to assign the thermodynamic driving force to a reaction.
"""
function _dg(model, rid_enzyme, rid_dg0, rid, θ; RT = 298.15 * 8.314e-3)
    rs = reaction_stoichiometry(model, rid)
    stoich = values(rs)
    mids = collect(keys(rs))
    midxs = Int.(indexin(mids, metabolites(model))) .+ length(rid_enzyme)
    if haskey(rid_dg0, rid)
        dg_val =
            rid_dg0[rid] + RT * sum(nu * log(θ[midx]) for (nu, midx) in zip(stoich, midxs))
        return 1.0 - exp(dg_val / RT)
    else # no kinetic info
        return 1.0
    end
end

"""
    $(TYPEDSIGNATURES)

Return a simplified version of `model` that contains only reactions (and the
associated metabolites) that are active, i.e. carry fluxes (from
`reaction_fluxes`) absolutely bigger than `atol`. All reactions are set
unidirectional based on `reaction_fluxes`. No gene information is copied.
"""
function prune_model(model::StandardModel, reaction_fluxes; atol = 1e-9, verbose = true)
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

    add_reactions!(pruned_model, rxns)
    add_metabolites!(pruned_model, mets)
    add_genes!(pruned_model, gs) #  no genes are added, unnecessary

    #: print some info about process
    if verbose
        rrn = n_reactions(model) - n_reactions(pruned_model)
        @info "Removed $rrn reactions."
    end

    return pruned_model
end
