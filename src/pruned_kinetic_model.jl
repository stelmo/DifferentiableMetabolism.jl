
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
function prune_fbc_model(model, ec_solution, flux_zero_tol, gene_zero_tol)

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

Build a pruned kinetic model using `model` as the base and adding kinetic
information using `reaction_isozymes`. Use `ec_solution` to prune away inactive
reactions, genes, metabolites, and isozymes. 

This model assumes enzyme kinetics are used, and thus
`gene_product_molar_masses` and `capacity` need to be supplied to limit the
supply of enzymes that can catalyze flux.

Note, a mixture of parameterized and non-parameterized structures can be
supplied. Finally, the model assumes that a turnover number like entity is
supplied in `reaction_isozymes`, but this can be any symbolic expression,
possibly including nonlinear kinetics. The only restriction is that these
kinetics may only depend on parameters.
"""
function build_pruned_kinetic_model(
    model,
    ec_solution,
    reaction_isozymes::Union{
        Dict{String,Dict{String,COBREXA.Isozyme}},
        Dict{String,Dict{String,ParameterIsozyme}},
    },
    gene_product_molar_masses::Union{Dict{String,Float64},Dict{String,Symbolics.Num}},
    capacity::Union{
        Vector{Tuple{String,Vector{String},Float64}},
        Float64,
        Symbolics.Num,
        Vector{Tuple{String,Vector{String},Symbolics.Num}},
    },
    flux_zero_tol,
    gene_zero_tol,
)
    pruned = prune_fbc_model(model, ec_solution, flux_zero_tol, gene_zero_tol)

    active_isozymes = Dict(
        k => iso for
        c in (ec_solution.isozyme_forward_amounts, ec_solution.isozyme_reverse_amounts)
        for (k, v) in c for (iso, val) in v if val > gene_zero_tol
    )

    pruned_reaction_isozymes = Dict(
        k => Dict(kk => vv for (kk, vv) in v if Symbol(kk) == active_isozymes[Symbol(k)]) for (k, v) in reaction_isozymes if haskey(active_isozymes, Symbol(k))
    )

    # build pruned model
    m = COBREXA.fbc_model_constraints(pruned)
    
    # less structure is required than with usual ec models, because the remaining reactiosn 
    m =
        m +
        :gene_product_amounts^ConstraintTrees.variables(
            keys = Symbol.(AbstractFBCModels.genes(pruned)),
            bounds = ConstraintTrees.Between(0, Inf),
        ) +
        :isozyme_amounts^COBREXA.isozyme_amount_variables(
            Symbol.(keys(pruned_reaction_isozymes)),
            rid -> Symbol.(keys(pruned_reaction_isozymes[string(rid)])),
        )

    # forward reactions, to correct sign of kcat
    fwd = [
        Symbol(rid) for rid in AbstractFBCModels.reactions(pruned) if
        pruned.reactions[rid].lower_bound >= 0
    ]

    kcat(rid, id) =
        getfield(
            pruned_reaction_isozymes[string(rid)][string(id)],
            rid in fwd ? :kcat_forward : :kcat_reverse,
        ) * (rid in fwd ? 1.0 : -1.0)

    m =
        m *
        :isozyme_flux_balance^ConstraintTrees.ConstraintTree(
            rid => ConstraintTrees.Constraint(
                sum(kcat(rid, id) * ri.value for (id, ri) in ris) - m.fluxes[rid].value,
                0.0,
            ) for (rid, ris) in m.isozyme_amounts
        ) *
        :gene_product_isozyme_balance^COBREXA.gene_product_isozyme_constraints(
            m.gene_product_amounts,
            tuple(m.isozyme_amounts),
            (rid, iid) -> COBREXA.maybemap(
                x -> [(Symbol(k), v) for (k, v) in x.gene_product_stoichiometry],
                COBREXA.maybeget(pruned_reaction_isozymes, string(rid), string(iid)),
            ),
        ) *
        :gene_product_capacity^(
            capacity isa Real ?
            ConstraintTrees.Constraint(
                value = sum(
                    gpa.value * gene_product_molar_masses[String(gp)] for
                    (gp, gpa) in m.gene_product_amounts
                ),
                bound = capacity isa Symbolics.Num ? ParameterBetween(0, capacity) :
                        ConstraintTrees.Between(0, capacity),
            ) :
            C.ConstraintTree(
                Symbol(id) => ConstraintTrees.Constraint(
                    value = sum(
                        m.gene_product_amounts[Symbol(gp)].value *
                        gene_product_molar_masses[gp] for gp in gps
                    ),
                    bound = limit isa Symbolics.Num ? ParameterBetween(0, limit) :
                            ConstraintTrees.Between(0, limit),
                ) for (id, gps, limit) in capacity_limits
            )
        )

    m
end

export build_pruned_kinetic_model
