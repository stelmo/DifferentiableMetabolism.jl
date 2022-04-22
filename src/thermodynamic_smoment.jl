"""
    $(TYPEDSIGNATURES)

Construct a [`DifferentiableModel`](@ref) from a [`COBREXA.SMomentModel`](@ref),
which includes kinetic as well as thermodynamic parameters. Add the standard
Gibbs free energy of reactions with `rid_dg0`, and the metabolite concentrations
that are to be used as parameters with `mid_concentration` (both arguments are
dictionaries mapping reaction or metabolite ids to values). See
`with_parameters(::SMomentModel,...)` for more information about restrictions to
the model types and arguments. 

Note, thermodynamic parameters require special attention. To ignore some
reactions when calculating the thermodynamic factor, include them in
`ignore_reaction_ids`. The units used for the Gibbs free energy are kJ/mol, and
concentrations should be in molar. Importantly, the standard Gibbs free energy
of reaction is assumed to be given *in the forward direction of the reaction* in
the model.
"""
function with_parameters(
    smm::SMomentModel,
    rid_enzyme::Dict{String,Enzyme},
    rid_dg0::Dict{String,Float64},
    mid_concentration::Dict{String,Float64};
    analytic_parameter_derivatives = x -> nothing,
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
)
    param_ids = [
        "k#" .* collect(keys(rid_enzyme))
        "c#" .* metabolites(smm)
    ]
    θ = [
        [x.kcat for x in values(rid_enzyme)]
        [mid_concentration[mid] for mid in metabolites(smm)]
    ]

    c, E, d, M, h, var_ids = _differentiable_thermodynamic_smoment_opt_problem(
        smm,
        rid_enzyme,
        rid_dg0;
        ϵ,
        atol,
        digits,
        RT,
        ignore_reaction_ids,
    )

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
    $(TYPEDSIGNATURES)

Return structures that will allow the most basic form of smoment to be solved.
No enzyme constraints allowed. Most effective enzyme is the only GRR. Assume
unidirectional reactions.
"""
function _differentiable_thermodynamic_smoment_opt_problem(
    smm::SMomentModel,
    rid_enzyme::Dict{String,Enzyme},
    rid_dg0::Dict{String,Float64};
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
)

    #: get irreverible stoichiometric matrix from model
    irrev_S = stoichiometry(smm.inner) * COBREXA._smoment_column_reactions(smm)

    #: make full rank
    S = DifferentiableMetabolism._remove_lin_dep_rows(irrev_S; ϵ, digits, atol)

    #: size of resultant model
    num_reactions = size(S, 2)
    num_eq_cons = size(S, 1)

    #: equality lhs
    E = _ -> S

    #: equality rhs
    d = _ -> spzeros(num_eq_cons) #TODO handle fixed variables

    #: objective inferred from model
    c = _ -> -objective(smm)

    #: coupling kcats and thermodynamics
    col_idxs, kcat_idxs =
        DifferentiableMetabolism._build_smoment_kcat_coupling(smm, rid_enzyme)

    kcat_thermo_coupling(θ) = sparsevec(
        col_idxs,
        [
            mw / (
                (θ[rid_idx]) * DifferentiableMetabolism._dg(
                    smm,
                    rid_enzyme,
                    rid_dg0,
                    rid,
                    θ;
                    RT,
                    ignore_reaction_ids,
                )
            ) for (rid, rid_idx, mw) in kcat_idxs
        ],
        n_reactions(smm),
    )

    #: inequality rhs
    Cp = coupling(smm.inner)
    M(θ) = [
        -1.0 * I(num_reactions)
        1.0 * I(num_reactions)
        Cp
        -kcat_thermo_coupling(θ)'
        kcat_thermo_coupling(θ)'
    ]

    #: inequality lhs
    xlb, xub = bounds(smm)
    clb, cub = coupling_bounds(smm.inner)
    h = _ -> [-xlb; xub; -clb; cub; 0; smm.total_enzyme_capacity]

    return c, E, d, M, h, reactions(smm)
end
