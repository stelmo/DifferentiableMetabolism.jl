"""
    with_parameters(
        gm::GeckoModel,
        rid_enzyme::Dict{String,Enzyme},
        rid_dg0::Dict{String,Float64},
        mid_concentration::Dict{String,Float64};
        scale_equality = false,
        analytic_parameter_derivatives = x -> nothing,
        ϵ = 1e-8,
        atol = 1e-12,
        digits = 8,
        RT = 298.15 * 8.314e-3,
        ignore_reaction_ids = [],
    )

Construct a [`DifferentiableModel`](@ref) from a [`COBREXA.GeckoModel`](@ref).
Each variable in `gm` is differentiated with respect to the kcats in the
dictionary `rid_enzyme`, which is a dictionary mapping reaction ids to
[`Enzyme`](@ref)s. Enzyme constraints are only taken with respect to the entries
of `rid_enzyme`. Incorporates thermodynamic constraints through `rid_dg0`. See
`with_parameters(::GeckoModel,...)` without the thermodynamic effects for
further details.

Note, thermodynamic parameters require special attention. To ignore some
reactions when calculating the thermodynamic factor, include them in
`ignore_reaction_ids`. The units used for the Gibbs free energy are kJ/mol, and
concentrations should be in molar. Importantly, the standard Gibbs free energy
of reaction is assumed to be given *in the forward direction of the reaction* in
the model.
"""
function with_parameters(
    gm::GeckoModel,
    rid_enzyme::Dict{String,Enzyme},
    rid_dg0::Dict{String,Float64},
    mid_concentration::Dict{String,Float64};
    scale_equality = false,
    analytic_parameter_derivatives = x -> nothing,
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

    c, _E, d, M, h, var_ids = _differentiable_thermodynamic_gecko_opt_problem(
        gm,
        rid_enzyme,
        rid_dg0;
        ϵ,
        atol,
        digits,
        RT,
        ignore_reaction_ids,
    )

    _make_differentiable_model(
        c,
        _E,
        d,
        M,
        h,
        θ,
        analytic_parameter_derivatives,
        param_ids,
        var_ids;
        scale_equality,
    )
end

"""
    _differentiable_thermodynamic_gecko_opt_problem(
        gm::GeckoModel,
        rid_enzyme::Dict{String,Enzyme},
        rid_dg0::Dict{String,Float64};
        ϵ = 1e-8,
        atol = 1e-12,
        digits = 8,
        RT = 298.15 * 8.314e-3,
        ignore_reaction_ids = [],
    )

Return optimization problem for a thermodynamic gecko problem, but in
differentiable format.
"""
function _differentiable_thermodynamic_gecko_opt_problem(
    gm::GeckoModel,
    rid_enzyme::Dict{String,Enzyme},
    rid_dg0::Dict{String,Float64};
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
                θ[ridx] * DifferentiableMetabolism._dg(
                    gm,
                    rid_enzyme,
                    rid_dg0,
                    rid,
                    θ;
                    RT,
                    ignore_reaction_ids,
                )
            ) for (rid, ridx, stoich) in kcat_rid_ridx_stoich
        ],
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
