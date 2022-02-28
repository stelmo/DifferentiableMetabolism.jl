"""
    remove_isozymes(model, reaction_kcats, protein_stoichiometry, gid_measurements)

Remove isozymes that are not expressed.
If multiple isozymes are expressed, pick one that has the highest expression.
"""
function remove_low_expressed_isozymes!(
    model::StandardModel,
    reaction_kcats,
    protein_stoichiometry,
    protein_masses,
    gid_measurements,
)

    for rid in reactions(model)
        if COBREXA._has_grr(model, rid)
            measured_proteins = Float64[]
            grrs = reaction_gene_association(model, rid)
            for (i, grr) in enumerate(grrs)

                push!(
                    measured_proteins,
                    sum(
                        map(
                            *,
                            protein_stoichiometry[rid][i],
                            [get(gid_measurements, gid, 0.0) for gid in grr],
                            [protein_masses[gid] for gid in grr],
                        ),
                    ),
                )
            end
            idx = argmax(measured_proteins)

            model.reactions[rid].grr = [grrs[idx]]
            reaction_kcats[rid] = [reaction_kcats[rid][idx]]
            protein_stoichiometry[rid] = [protein_stoichiometry[rid][idx]]
        end
    end

    curated_gids = String[]
    for rid in reactions(model)
        if COBREXA._has_grr(model, rid)
            for grr in reaction_gene_association(model, rid)
                append!(curated_gids, grr)
            end
        end
    end
    rm_gids = setdiff(genes(model), curated_gids)
    delete!(model.genes, rm_gids) # remove genes that were deleted
    return nothing
end

"""
    prune_model(
        model::StandardModel,
        reaction_fluxes;
        rtol = 1e-9,
    )

Return a simplified version of `model` that contains only reactions 
(and the associated metabolites and genes) that are active, i.e. carry 
fluxes (from `reaction_fluxes`) absolutely bigger than `rtol`. 

Assume:
1) Model does not have isozymes, isozyme with the fastest kcat has been preselected
"""
function prune_model(model::StandardModel, reaction_fluxes; rtol = 1e-9)
    pruned_model = StandardModel("pruned_model")

    rxns = Vector{Reaction}()
    mets = Vector{Metabolite}()
    gs = Vector{Gene}()

    for rid in reactions(model)
        isapprox(abs(reaction_fluxes[rid]), 0; atol = rtol) && continue

        rxn = deepcopy(model.reactions[rid])
        if reaction_fluxes[rid] > 0
            rxn.lb = max(0, rxn.lb)
        else
            rxn.ub = min(0, rxn.ub)
        end
        push!(rxns, rxn)

        rs = reaction_stoichiometry(model, rid)
        for mid in keys(rs)
            push!(mets, model.metabolites[mid])
        end

        grrs = reaction_gene_association(model, rid)
        isnothing(grrs) && continue
        for grr in grrs
            for gid in grr
                push!(gs, model.genes[gid])
            end
        end
    end

    add_reactions!(pruned_model, rxns)
    add_metabolites!(pruned_model, mets)
    add_genes!(pruned_model, gs)

    #: print some info about process
    rrn = n_reactions(model) - n_reactions(pruned_model)
    @info "Removed $rrn reactions."
    @info "Pruned model has $(n_reactions(pruned_model)) reactions, $(n_metabolites(pruned_model)) metabolites, and $(n_genes(pruned_model)) genes."

    return pruned_model
end

"""
    _remove_lin_dep_rows(A; ϵ=1e-8)

Remove linearly dependent rows using reduced row echelon form. Adjust sensitivity to 
numerical issues with `ϵ`.
"""
function _remove_lin_dep_rows(A; ϵ = 1e-8)
    #TODO this method is suboptimal and can be improved for numerical stability
    #TODO improve RowEchelon, SVD does not work due to column reordering
    rA = rref!(copy(Array(A)), ϵ)
    idxs = Int[]
    for i = 1:size(rA, 1)
        if !all(rA[i, :] .== 0) # remove rows of all zero
            push!(idxs, i)
        end
    end
    remrows = size(A, 1) - length(idxs)
    @info "Removed $remrows rows!"
    return rA[idxs, :]
end

function in_another_grr(model, current_rid, current_gid)
    for rid in reactions(model)
        current_rid == rid && continue 
        !COBREXA._has_grr(model, rid) && continue 
        for grr in reaction_gene_association(model, rid)
            for gid in grr 
                current_gid == gid && return true
            end
        end
    end
    false
end