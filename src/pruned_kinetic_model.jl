function prune_fbc_model(model, ec_solution, zerotol)

    rids = [string(k) for (k, v) in ec_solution.fluxes if abs(v) > zerotol]
    gids = [string(k) for (k, v) in ec_solution.gene_product_amounts if abs(v) > zerotol]
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
    zerotol,
)
    pruned = prune_fbc_model(model, ec_solution, zerotol)

    active_isozymes = Dict(
        k => iso for
        c in (ec_solution.isozyme_forward_amounts, ec_solution.isozyme_reverse_amounts)
        for (k, v) in c for (iso, val) in v if val > zerotol
    )

    pruned_reaction_isozymes = Dict(
        k => Dict(kk => vv for (kk, vv) in v if Symbol(kk) == active_isozymes[Symbol(k)]) for
        (k, v) in reaction_isozymes if haskey(active_isozymes, Symbol(k))
    )

    # build pruned model
    m = COBREXA.fbc_model_constraints(pruned)

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

    # connect all parts with constraints
    fwd = [Symbol(rid) for rid in AbstractFBCModels.reactions(pruned) if pruned.reactions[rid].lower_bound >= 0]
    
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
                    (gp, gpa) in constraints.gene_product_amounts
                ),
                bound = capacity isa Symbolics.Num ? ParameterBetween(0, capacity) :
                        ConstraintTrees.Between(0, capacity),
            ) :
            C.ConstraintTree(
                Symbol(id) => ConstraintTrees.Constraint(
                    value = sum(
                        constraints.gene_product_amounts[Symbol(gp)].value *
                        gene_product_molar_masses[gp] for gp in gps
                    ),
                    bound = limit isa Symbolics.Num ? ParameterBetween(0, limit) :
                            ConstraintTrees.Between(0, limit),
                ) for (id, gps, limit) in capacity_limits
            )
        )

    m
end

export build_pruned_model
