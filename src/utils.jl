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

# """
#     $(TYPEDSIGNATURES)

# Return a simplified version of `model` that contains only reactions 
# (and the associated metabolites and genes) that are active, i.e. carry 
# fluxes (from `reaction_fluxes`) absolutely bigger than `rtol`. All 
# reactions are unidirectional. 
# """
# function prune_model(
#     model::StandardModel, 
#     reaction_fluxes; 
#     rtol = 1e-9
# )
#     pruned_model = StandardModel("pruned_model")

#     rxns = Vector{Reaction}()
#     mets = Vector{Metabolite}()
#     gs = Vector{Gene}()

#     for rid in reactions(model)
#         isapprox(abs(reaction_fluxes[rid]), 0; atol = rtol) && continue

#         rxn = deepcopy(model.reactions[rid])
#         if reaction_fluxes[rid] > 0
#             rxn.lb = max(0, rxn.lb)
#         else
#             rxn.ub = min(0, rxn.ub)
#         end
#         push!(rxns, rxn)

#         rs = reaction_stoichiometry(model, rid)
#         for mid in keys(rs)
#             push!(mets, model.metabolites[mid])
#         end

#         !has_reaction_grr(model, rid) && continue
#         for grr in reaction_gene_association(model, rid)
#             for gid in grr
#                 push!(gs, model.genes[gid])
#             end
#         end
#     end

#     add_reactions!(pruned_model, rxns)
#     add_metabolites!(pruned_model, mets)
#     add_genes!(pruned_model, gs)

#     #: print some info about process
#     rrn = n_reactions(model) - n_reactions(pruned_model)
#     @info "Removed $rrn reactions."

#     return pruned_model
# end

# function in_another_grr(model, current_rid, current_gid)
#     for rid in reactions(model)
#         current_rid == rid && continue
#         !COBREXA._has_grr(model, rid) && continue
#         for grr in reaction_gene_association(model, rid)
#             for gid in grr
#                 current_gid == gid && return true
#             end
#         end
#     end
#     false
# end



# """
# """
# function qp_objective_measured(
#     rids,
#     gids,
#     obs_v_dict,
#     obs_e_dict;
#     vtol = 1e-3,
#     etol = 1e-3,
#     reg = 1e-1,
# )
#     n_vars = length(gids) + length(rids)
#     c = zeros(n_vars)
#     q = zeros(n_vars)
#     n = 0

#     for (i, rid) in enumerate(rids)
#         if !haskey(obs_v_dict, rid) || abs(obs_v_dict[rid]) < vtol
#             q[i] = reg
#         else
#             c[i] = -1.0 / abs(obs_v_dict[rid]) #! fluxes are positive in model
#             q[i] = 1.0 / obs_v_dict[rid]^2
#             n += 1
#         end
#     end
#     k = length(rids)
#     for (i, gid) in enumerate(gids)
#         if !haskey(obs_e_dict, gid) || abs(obs_e_dict[gid]) < etol
#             q[k+i] = reg
#         else
#             c[k+i] = -1.0 / obs_e_dict[gid]
#             q[k+i] = 1.0 / obs_e_dict[gid]^2
#             n += 1
#         end
#     end

#     return spdiagm(q), sparse(c), n
# end
