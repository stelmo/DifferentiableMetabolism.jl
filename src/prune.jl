
#=
Copyright (c) 2025, Heinrich-Heine University Duesseldorf
Copyright (c) 2025, University of Luxembourg

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

"""
$(TYPEDSIGNATURES)

Return a freshly allocated dictionary mapping reaction IDs to isozyme structs.
Note, all reactions running in reverse have their kcats swapped, since
[`prune_model`](@ref) changes all fluxes to be positive.

# Example
```
prune_reaction_isozymes(reaction_isozymes, solution.isozyme_forward_amounts, solution.isozyme_reverse_amounts, solution.fluxes, 1e-6)
```
"""
function prune_reaction_isozymes(
    reaction_isozymes,
    solution_isozyme_forward_amounts,
    solution_isozyme_reverse_amounts,
    solution_fluxes,
    gene_zero_tol::Float64,
)
    active_isozymes = [
        rid => iso_id for iso_amounts in
        (solution_isozyme_forward_amounts, solution_isozyme_reverse_amounts) for
        (rid, iso_vals) in iso_amounts for
        (iso_id, iso_val) in iso_vals if iso_val > gene_zero_tol
    ]

    @assert length(x.first for x in active_isozymes) ==
            length(unique(x.first for x in active_isozymes))

    remove_unused_kcat(isostruct, flux) = begin
        x = deepcopy(isostruct)
        if flux < 0
            x.kcat_forward = x.kcat_reverse # swap to make all reactions forward, also in pruned model
        end
        x.kcat_reverse = nothing

        x
    end

    Dict(  # pruned_reaction_isozymes
        string(rid) => Dict(
            string(iso_id) => remove_unused_kcat(
                reaction_isozymes[string(rid)][string(iso_id)],
                solution_fluxes[rid],
            ),
        ) for (rid, iso_id) in active_isozymes
    )
end

"""
$(TYPEDSIGNATURES)

Prune away reactions, metabolites, and genes from a `model` using `solution`.
Additionally, adjust `reaction_isozymes` using [`prune_reaction_isozymes`](@ref)
to account for the changed directions of the fluxes.

Fluxes and gene product concentrations smaller than `flux_zero_tol`,
`gene_zero_tol` are removed. Metabolites that do not take part in the remaining
reactions are also removed.

Note, in COBREXA reverse variables are always derived quantities, hence making
everything forward orientated ensures that the variables of the whole model are
just these.
"""
function prune_model(
    model::A.CanonicalModel.Model,
    solution_fluxes,
    solution_gene_product_amounts,
    reaction_isozymes,
    solution_isozyme_forward_amounts,
    solution_isozyme_reverse_amounts,
    flux_zero_tol::Float64,
    gene_zero_tol::Float64,
)

    rids = [string(k) for (k, v) in solution_fluxes if abs(v) > flux_zero_tol]
    gids = [string(k) for (k, v) in solution_gene_product_amounts if abs(v) > gene_zero_tol]
    mids = Set(mid for rid in rids for mid in keys(A.reaction_stoichiometry(model, rid)))

    d_rids = setdiff(A.reactions(model), rids)
    d_mids = setdiff(A.metabolites(model), mids)
    d_gids = setdiff(A.genes(model), gids)

    pruned = convert(A.CanonicalModel.Model, model)

    for rid in d_rids
        delete!(pruned.reactions, rid)
    end
    for mid in d_mids
        delete!(pruned.metabolites, mid)
    end
    for gid in d_gids
        delete!(pruned.genes, gid)
    end

    # set bounds - should all be unidirectional now
    for rid in [rid for rid in rids if solution_fluxes[Symbol(rid)] > 0]
        pruned.reactions[rid].lower_bound = max(pruned.reactions[rid].lower_bound, 0)
    end
    # additionally, make all reverse flux reactions forward flux - NB: assumes the kcats have been swapped as well - simplifies things downstream
    for rid in [rid for rid in rids if solution_fluxes[Symbol(rid)] < 0]
        pruned.reactions[rid].stoichiometry =
            Dict(k => -1 * v for (k, v) in pruned.reactions[rid].stoichiometry)
        lb = pruned.reactions[rid].lower_bound
        ub = pruned.reactions[rid].upper_bound
        pruned.reactions[rid].lower_bound = ub > 0 ? min(0, abs(ub)) : max(0, abs(ub))
        pruned.reactions[rid].upper_bound = abs(lb)
    end

    pruned,
    prune_reaction_isozymes(
        reaction_isozymes,
        solution_isozyme_forward_amounts,
        solution_isozyme_reverse_amounts,
        solution_fluxes,
        gene_zero_tol,
    )
end

export prune_model
