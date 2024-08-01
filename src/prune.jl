
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

Changes from copied code are indicated.
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
        mid in keys(reaction_stoichiometry(model, rid))
    )

    d_rids = setdiff(reactions(model), rids)
    d_mids = setdiff(metabolites(model), mids)
    d_gids = setdiff(genes(model), gids)

    pruned = convert(AC.Model, model)

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

function prune_reaction_isozymes(reaction_isozymes, ec_solution, flux_zero_tol, gene_zero_tol)
    active_isozymes = Dict(
        k => iso for c in (ec_solution.isozyme_forward_amounts, ec_solution.isozyme_reverse_amounts) for (k, v) in c for (iso, val) in v if val > gene_zero_tol
    )

    Dict(  # pruned_reaction_isozymes
        k => Dict(kk => vv for (kk, vv) in v if Symbol(kk) == active_isozymes[Symbol(k)]) for (k, v) in reaction_isozymes if haskey(active_isozymes, Symbol(k))
    )
end

export prune_reaction_isozymes

function prune_gene_product_molar_masses(gene_product_molar_masses, ec_solution, flux_zero_tol, gene_zero_tol)
    Dict(
        string(k) => gene_product_molar_masses[string(k)] for (k, v) in ec_solution.gene_product_amounts if v >= gene_zero_tol
    )
end

export prune_gene_product_molar_masses
