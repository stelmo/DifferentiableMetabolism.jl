"""
    with_parameters(
        smm::SMomentModel,
        rid_enzyme::Dict{String,Enzyme};
        analytic_parameter_derivatives = x -> nothing,
        ϵ = 1e-8,
        atol = 1e-12,
        digits = 8,
    )

Construct a [`DifferentiableModel`](@ref) from a [`COBREXA.SMomentModel`](@ref).
Each variable in `smm` is differentiated with respect to the kcats in the
dictionary `rid_enzyme`, which is a dictionary mapping reaction ids to
[`Enzyme`](@ref)s.Enzyme constraints are only taken with respect to the entries
of `rid_enzyme`.

The analytic derivative of the optimality conditions with respect to the
parameters can be supplied through `analytic_parameter_derivatives`. Internally,
`ϵ`, `atol`, and `digits`, are forwarded to [`_remove_lin_dep_rows`](@ref). 

Note, to ensure differentiability, preprocessing of the model is required. In short,
only an active solution may be differentiated, this required that:
- the model does not possess any isozymes
- all the reactions should be unidirectinal
- the kcats in `rid_enzyme` are for the appropriate direction used in the model
- all rids in `rid_enzyme` are used in the model
"""
function with_parameters(
    smm::SMomentModel,
    rid_enzyme::Dict{String,Enzyme};
    analytic_parameter_derivatives = x -> nothing,
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
)
    param_ids = "k#" .* collect(keys(rid_enzyme))
    θ = [x.kcat for x in values(rid_enzyme)]

    c, E, d, M, h, var_ids =
        _differentiable_smoment_opt_problem(smm, rid_enzyme; ϵ, atol, digits)

    return DifferentiableModel(
        _ -> spzeros(length(var_ids), length(var_ids)),
        c,
        E,
        d,
        M,
        h,
        θ,
        analytic_parameter_derivatives,
        param_ids,
        var_ids,
    )
end

"""
    _differentiable_smoment_opt_problem(
        smm::SMomentModel,
        rid_enzyme::Dict{String,Enzyme};
        ϵ = 1e-8,
        atol = 1e-12,
        digits = 8,
    )

Return structures that will allow the most basic form of smoment to be solved.
No enzyme constraints allowed. Assume preprocessing changes model such that most
effective enzyme is the only GRR.
"""
function _differentiable_smoment_opt_problem(
    smm::SMomentModel,
    rid_enzyme::Dict{String,Enzyme};
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
)

    #: get irreverible stoichiometric matrix from model
    irrev_S = stoichiometry(smm.inner) * COBREXA._smoment_column_reactions(smm)

    #: make full rank
    S = DifferentiableMetabolism._remove_lin_dep_rows(irrev_S; ϵ, digits, atol)

    #: size of resultant model
    num_reactions = size(S, 2)
    num_metabolites = size(S, 1)

    #: equality lhs
    E = _ -> S

    #: equality rhs
    d = _ -> spzeros(num_metabolites) #TODO handle fixed variables

    #: objective inferred from model
    c = _ -> -objective(smm)

    #: coupling to kcats
    col_idxs, kcat_idxs = _build_smoment_kcat_coupling(smm, rid_enzyme)

    kcat_coupling(θ) =
        sparsevec(col_idxs, [mw / θ[idx] for (_, idx, mw) in kcat_idxs], num_reactions)

    #: inequality rhs
    Cp = coupling(smm.inner)
    M(θ) = [
        -1.0 * I(num_reactions)
        1.0 * I(num_reactions)
        Cp
        -kcat_coupling(θ)'
        kcat_coupling(θ)'
    ]

    #: inequality lhs
    xlb, xub = bounds(smm)
    clb, cub = coupling_bounds(smm.inner)
    h = _ -> [-xlb; xub; -clb; cub; 0; smm.total_enzyme_capacity]

    return c, E, d, M, h, reactions(smm)
end

"""
    _build_smoment_kcat_coupling(smm::SMomentModel, rid_enzyme)

Helper function to build kcat coupling in parametric form for smoment problems.
"""
function _build_smoment_kcat_coupling(smm::SMomentModel, rid_enzyme)
    kcat_original_rid_order = String[]
    col_idxs = Int[]
    mws = Float64[]
    for (col_idx, rid) in enumerate(reactions(smm))
        original_rid = string(first(split(rid, "#")))

        # skip these entries
        !haskey(rid_enzyme, original_rid) && continue

        # these entries have kcats, only one GRR by assumption
        mw = sum([
            pstoich * rid_enzyme[original_rid].gene_product_mass[gid] for
            (gid, pstoich) in rid_enzyme[original_rid].gene_product_count
        ])
        push!(kcat_original_rid_order, original_rid)
        push!(col_idxs, col_idx)
        push!(mws, mw)
    end

    kcat_idxs = [
        (rid, first(indexin([rid], collect(keys(rid_enzyme)))), mw) for
        (rid, mw) in zip(kcat_original_rid_order, mws)
    ]

    return col_idxs, kcat_idxs
end
