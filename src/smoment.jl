"""
    $(TYPEDSIGNATURES)

Construct a [`DifferentiableModel`](@ref) from a [`SMomentModel`]. Each variable
in `gm` is differentiated with respect to the kcats in the dictionary
`rid_enzyme`, which is a dictionary mapping reaction ids to [`Enzyme`](@ref)s.
Only an active solution may be differentiated, i.e. the model should not possess
any isozymes, all reactions are unidirectinal, and the kcats in `rid_enzyme`
should be for the appropriate direction used in the model. Use
[`prune_model`](@ref) to ensure that an appropriate model is differentiated.  
"""
function with_parameters(
    smm::SMomentModel, 
    rid_enzyme::Dict{String, Enzyme};
    analytic_parameter_derivatives= x -> nothing,
    ϵ = 1e-8,
    stoich_digits_round = 8
)   
    param_ids = "k#" .* collect(keys(rid_enzyme))
    θ = [x.kcat for x in values(rid_enzyme)]

    c, Sf, d, M, hf, var_ids = differentiable_smoment_opt_problem(
        smm,
        rid_enzyme;
        ϵ,
        stoich_digits_round,
    )

    return DifferentiableModel(
        spzeros(0,0),
        c,
        Sf,
        d,
        M, 
        hf,
        θ,
        analytic_parameter_derivatives,
        param_ids,
        var_ids,
    )
end

"""
    $(TYPEDSIGNATURES)

Return structures that will allow the most basic form of smoment to be solved.
No enzyme constraints allowed. Assume preprocessing changes model such that most
effective enzyme is the only GRR.
"""
function differentiable_smoment_opt_problem(
    model::SMomentModel,
    rid_enzyme::Dict{String, Enzyme};
    ϵ = 1e-8,
    stoich_digits_round = 8,
)

    #: get irreverible stoichiometric matrix from model
    irrev_S = model.S[1:length(model.metabolites), 1:length(model.irrev_reaction_ids)]

    #: make full rank
    S = round.(_remove_lin_dep_rows(irrev_S; ϵ), digits = stoich_digits_round)

    #: size of resultant model
    num_reactions = size(S, 2)
    num_metabolites = size(S, 1)
    num_vars = num_reactions + 1

    #: equality lhs
    kcat_original_rid_order = String[]
    col_idxs = Int[]
    neg_mws = Float64[]
    for (col_idx, rid) in enumerate(model.irrev_reaction_ids)
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        !haskey(rid_enzyme, original_rid) && continue

        # these entries have kcats, only one GRR by assumption
        mw = sum([pstoich * rid_enzyme[original_rid].molar_masses[gid] for (gid, pstoich) in rid_enzyme[original_rid].stoichiometry])
        push!(kcat_original_rid_order, original_rid)
        push!(col_idxs, col_idx)
        push!(neg_mws, -mw)
    end

    kcat_idxs = [
        (first(indexin([rid], collect(keys(rid_enzyme)))), nmw) for
        (rid, nmw) in zip(kcat_original_rid_order, neg_mws)
    ]

    fSe(θ) = sparsevec(col_idxs, [nmw / θ[x] for (x, nmw) in kcat_idxs], num_reactions)

    Ef(θ) = [
        Array(S) zeros(num_metabolites, 1)
        fSe(θ)' 1.0
    ]

    #: equality rhs
    d = spzeros(num_metabolites + 1)

    #: objective inferred from model
    c = -objective(model)

    #: inequality constraints
    xlb, xub = bounds(model)
    M = [
        -1.0 * I(num_vars)
        1.0 * I(num_vars)
    ]

    h = [-xlb; xub]

    hf(x) = h

    return c, Ef, d, M, hf, model.irrev_reaction_ids
end
