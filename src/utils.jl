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

function rescale(mat, vec; verbose = false, atol = 1e-8)
    max_coeff_range, _ = check_scaling(mat; atol)
    verbose && println("Coefficient range: ", max_coeff_range)
    lb = -round(max_coeff_range, RoundUp) / 2.0
    ub = round(max_coeff_range, RoundUp) / 2.0

    rsmat = similar(mat)
    rsvec = similar(vec)
    for (j, row) in enumerate(eachrow(mat))
        llv, luv = extrema(log10 ∘ abs, row)
        if lb <= llv && luv <= ub
            sf = 0.0
        else # scale to upper bound
            sf = ub - luv
        end
        rsmat[j, :] .= row .* 10^sf
        rsvec[j] = vec[j] * 10^sf
    end

    return rsmat, rsvec
end

desc(x; atol = 1e-8) = abs(x) > atol && log10(abs(x))

function check_scaling(mat; atol = 1e-8)
    best_case = maximum(
        maximum(desc.(mat; atol), dims = 2)[:] - minimum(desc.(mat; atol), dims = 2)[:],
    )
    rlb, rub = log10.(extrema(filter(x -> x > atol, abs.(mat))))
    worst_case = rub - rlb
    return best_case, worst_case
end

