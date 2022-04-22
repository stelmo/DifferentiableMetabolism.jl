"""
    with_parameters(
        gm::GeckoModel,
        rid_enzyme::Dict{String,Enzyme};
        analytic_parameter_derivatives = x -> nothing,
        ϵ = 1e-8,
        atol = 1e-12,
        digits = 8,
    )

Construct a [`DifferentiableModel`](@ref) from a [`COBREXA.GeckoModel`](@ref).
Each variable in `gm` is differentiated with respect to the kcats in the
dictionary `rid_enzyme`, which is a dictionary mapping reaction ids to
[`Enzyme`](@ref)s. Enzyme constraints are only taken with respect to the entries
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
    gm::GeckoModel,
    rid_enzyme::Dict{String,Enzyme};
    analytic_parameter_derivatives = x -> nothing,
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
)
    param_ids = "k#" .* collect(keys(rid_enzyme))
    θ = [x.kcat for x in values(rid_enzyme)]

    c, E, d, M, h, var_ids =
        _differentiable_gecko_opt_problem(gm, rid_enzyme; ϵ, atol, digits)

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
    _differentiable_gecko_opt_problem(
        gm::GeckoModel,
        rid_enzyme::Dict{String,Enzyme};
        ϵ = 1e-8,
        atol = 1e-12,
        digits = 8,
    )

Return optimization problem for gecko problem, but in differentiable format.
"""
function _differentiable_gecko_opt_problem(
    gm::GeckoModel,
    rid_enzyme::Dict{String,Enzyme};
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
)
    #: get irreverible stoichiometric matrix from model
    irrev_S = stoichiometry(gm.inner) * COBREXA._gecko_reaction_column_reactions(gm)

    #: make full rank
    S = DifferentiableMetabolism._remove_lin_dep_rows(irrev_S; ϵ, digits, atol)

    #: size of resultant model
    num_reactions = size(S, 2)
    num_genes = n_genes(gm)
    num_metabolites = size(S, 1)
    num_vars = num_reactions + num_genes

    #: equality lhs
    E_components, kcat_rid_ridx_stoich =
        DifferentiableMetabolism._build_gecko_equality_enzyme_constraints(gm, rid_enzyme)

    Se(θ) = sparse(
        E_components.row_idxs,
        E_components.col_idxs,
        [stoich / θ[idx] for (_, idx, stoich) in kcat_rid_ridx_stoich],
        num_genes,
        num_reactions,
    )

    E(θ) = [
        S spzeros(num_metabolites, num_genes)
        Se(θ) I(num_genes)
    ]

    #: equality rhs
    d = _ -> spzeros(num_metabolites + num_genes) #TODO handle fixed variables

    #: objective inferred from model
    c = _ -> -objective(gm)

    #: inequality constraints
    Cp = coupling(gm)
    clb, cub = coupling_bounds(gm)
    xlb, xub = bounds(gm)
    M = _ -> [
        -I(num_vars)
        I(num_vars)
        -Cp
        Cp
    ]

    h = _ -> [-xlb; xub; -clb; cub]

    return c, E, d, M, h, [reactions(gm); genes(gm)]
end

"""
    _add_gecko_enzyme_variable_as_function(
        rid_enzyme,
        original_rid,
        E_components,
        col_idx,
        gene_ids,
    )

Helper function to add an column into the enzyme stoichiometric matrix
parametrically.
"""
function _add_gecko_enzyme_variable_as_function(
    rid_enzyme,
    original_rid,
    E_components,
    col_idx,
    gene_ids,
)
    for (gid, stoich) in rid_enzyme[original_rid].gene_product_count
        push!(E_components.row_idxs, first(indexin([gid], gene_ids)))
        push!(E_components.col_idxs, col_idx)
        push!(E_components.coeff_tuple, (-stoich, original_rid))
    end
end

"""
    _build_gecko_equality_enzyme_constraints(gm::GeckoModel, rid_enzyme)

Helper function to build the equality enzyme constraints.
"""
function _build_gecko_equality_enzyme_constraints(gm::GeckoModel, rid_enzyme)
    E_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeff_tuple = Vector{Tuple{Float64,String}}(),
    )

    gids = genes(gm)
    for (col_idx, rid) in enumerate(reactions(gm))
        original_rid = string(first(split(rid, "#")))

        !haskey(rid_enzyme, original_rid) && continue

        #: add all entries to column of matrix
        DifferentiableMetabolism._add_gecko_enzyme_variable_as_function(
            rid_enzyme,
            original_rid,
            E_components,
            col_idx,
            gids,
        )
    end

    kcat_rid_ridx_stoich = [
        (rid, first(indexin([rid], collect(keys(rid_enzyme)))), stoich) for
        (stoich, rid) in E_components.coeff_tuple
    ]

    return E_components, kcat_rid_ridx_stoich
end
