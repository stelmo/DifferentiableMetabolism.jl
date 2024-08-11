
#=
Copyright (c) 2023, Heinrich-Heine University Duesseldorf
Copyright (c) 2023, University of Luxembourg

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

Prune away reactions, metabolites, and genes from a `model` using `ec_solution`,
which is the result of an enzyme constrained kinetic simulation. Fluxes and gene
product concentrations smaller than `flux_zero_tol`, `gene_zero_tol` are
removed. Metabolites that do not take part in the remaining reactions are also
removed.
"""
function prune_model(model, ec_solution, flux_zero_tol, gene_zero_tol)

    rids = [string(k) for (k, v) in ec_solution.fluxes if abs(v) > flux_zero_tol]
    gids =
        [string(k) for (k, v) in ec_solution.gene_product_amounts if abs(v) > gene_zero_tol]
    mids = Set(
        mid for rid in rids for
        mid in keys(AbstractFBCModels.reaction_stoichiometry(model, rid))
    )

    d_rids = setdiff(AbstractFBCModels.reactions(model), rids)
    d_mids = setdiff(AbstractFBCModels.metabolites(model), mids)
    d_gids = setdiff(AbstractFBCModels.genes(model), gids)

    pruned = convert(AbstractFBCModels.CanonicalModel.Model, model)

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
    for rid in [rid for rid in rids if ec_solution.fluxes[Symbol(rid)] > 0]
        pruned.reactions[rid].lower_bound = max(pruned.reactions[rid].lower_bound, 0)
    end
    for rid in [rid for rid in rids if ec_solution.fluxes[Symbol(rid)] < 0]
        pruned.reactions[rid].upper_bound = min(0, pruned.reactions[rid].upper_bound)
    end

    pruned
end

export prune_model

"""
$(TYPEDSIGNATURES)

TODO
"""
function prune_reaction_isozymes(reaction_isozymes, ec_solution, gene_zero_tol)
    active_isozymes = [
        rid => iso_id for iso_amounts in
        (ec_solution.isozyme_forward_amounts, ec_solution.isozyme_reverse_amounts) for
        (rid, iso_vals) in iso_amounts for
        (iso_id, iso_val) in iso_vals if iso_val > gene_zero_tol
    ]

    @assert length(x.first for x in active_isozymes) ==
            length(unique(x.first for x in active_isozymes))

    remove_unused_kcat(isostruct, flux) = begin
        x = deepcopy(isostruct)
        if flux < 0
            x.kcat_forward = nothing
        else
            x.kcat_reverse = nothing
        end
        x
    end

    Dict(  # pruned_reaction_isozymes
        string(rid) => Dict(
            string(iso_id) => remove_unused_kcat(
                reaction_isozymes[string(rid)][string(iso_id)],
                ec_solution.fluxes[rid],
            ),
        ) for (rid, iso_id) in active_isozymes
    )
end

export prune_reaction_isozymes

"""
$(TYPEDSIGNATURES)

TODO
"""
function prune_gene_product_molar_masses(
    gene_product_molar_masses,
    ec_solution,
    gene_zero_tol,
)
    Dict(
        string(k) => gene_product_molar_masses[string(k)] for
        (k, v) in ec_solution.gene_product_amounts if v >= gene_zero_tol
    )
end

export prune_gene_product_molar_masses
