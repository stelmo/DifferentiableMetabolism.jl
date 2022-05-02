"""
$(TYPEDSIGNATURES)

Construct a [`DifferentiableModel`](@ref) from a [`COBREXA.GeckoModel`](@ref).
Each variable in `gm` is differentiated with respect to the kcats in the
dictionary `rid_enzyme`, which is a dictionary mapping reaction ids to
[`Enzyme`](@ref)s. Enzyme constraints are only taken with respect to the entries
of `rid_enzyme`. Incorporates thermodynamic and saturation constraints through
`rid_dg0` and `rid_km`. See `with_parameters(::GeckoModel,...)` without 
the thermodynamic and saturation effects for further details.
"""
function with_parameters(
    gm::GeckoModel,
    rid_enzyme::Dict{String,Enzyme},
    rid_dg0::Dict{String,Float64},
    rid_km::Dict{String,Dict{String,Float64}},
    mid_concentration::Dict{String,Float64};
    scale_equality = false,
    scale_inequality = false,
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
)
    param_ids = [
        "k#" .* collect(keys(rid_enzyme))
        "c#" .* metabolites(gm.inner)
    ]
    θ = [
        [x.kcat for x in values(rid_enzyme)]
        [mid_concentration[mid] for mid in metabolites(gm.inner)]
    ]

    c, E, d, M, h, var_ids = _differentiable_michaelis_menten_gecko_opt_problem(
        gm,
        rid_enzyme,
        rid_dg0,
        rid_km;
        ϵ,
        atol,
        digits,
        RT,
        ignore_reaction_ids,
    )

    _make_differentiable_model(
        c,
        E,
        d,
        M,
        h,
        θ,
        var_ids,
        param_ids;
        scale_equality,
        scale_inequality,
    )
end

"""
$(TYPEDSIGNATURES)

Return optimization problem where thermodynamic and saturation effects are
incorporated into the gecko problem, but in differentiable format.
"""
function _differentiable_michaelis_menten_gecko_opt_problem(
    gm::GeckoModel,
    rid_enzyme::Dict{String,Enzyme},
    rid_dg0::Dict{String,Float64},
    rid_km::Dict{String,Dict{String,Float64}};
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
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
        [
            stoich / (
                θ[ridx] *
                DifferentiableMetabolism._dg(
                    gm,
                    rid_enzyme,
                    rid_dg0,
                    rid,
                    θ;
                    RT,
                    ignore_reaction_ids,
                ) *
                DifferentiableMetabolism._saturation(
                    gm,
                    rid_enzyme,
                    rid_km,
                    rid,
                    mangled_rid,
                    θ;
                    ignore_reaction_ids,
                )
            ) for (mangled_rid, (rid, ridx, stoich)) in
            zip(reactions(gm)[E_components.col_idxs], kcat_rid_ridx_stoich)
        ],
        num_genes,
        num_reactions,
    )

    E(θ) = [
        S spzeros(num_metabolites, num_genes)
        Se(θ) spdiagm(fill(1.0, num_genes))
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
        -spdiagm(fill(1.0, num_vars))
        spdiagm(fill(1.0, num_vars))
        -Cp
        Cp
    ]

    h = _ -> [-xlb; xub; -clb; cub]

    return c, E, d, M, h, [reactions(gm); genes(gm)]
end
