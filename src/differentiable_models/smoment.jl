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
    rid_dg0 = Dict{String,Float64}(),
    mid_concentration = Dict{String,Float64}(),
    scale_equality = false,
    scale_inequality = false,
    ϵ = 1e-8,
    atol = 1e-12,
    digits = 8,
    RT = 298.15 * 8.314e-3,
    ignore_reaction_ids = [],
    ignore_metabolite_ids = [],
)
    if isempty(rid_dg0)
        # no concentration parameters
        param_ids = "k#" .* collect(keys(rid_enzyme))
        θ = [x.kcat for x in values(rid_enzyme)]
    else
        # has concentration parameters
        isempty(mid_concentration) && throw(error("Missing metabolite concentrations."))
        param_ids = [
            "k#" .* collect(keys(rid_enzyme))
            "c#" .* metabolites(smm)
        ]
        θ = [
            [x.kcat for x in values(rid_enzyme)]
            [mid_concentration[mid] for mid in metabolites(smm)]
        ]
    end

    c, E, d, M, h, var_ids = _differentiable_thermodynamic_smoment_opt_problem(
        smm,
        rid_enzyme,
        rid_dg0;
        ϵ,
        atol,
        digits,
        RT,
        ignore_reaction_ids,
        ignore_metabolite_ids,
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
    ignore_metabolite_ids = [],
)

    #: get irreverible stoichiometric matrix from model
    irrev_S = stoichiometry(smm.inner) * COBREXA._smoment_column_reactions(smm)

    #: make full rank
    S = sparse(DifferentiableMetabolism._remove_lin_dep_rows(irrev_S; ϵ, digits, atol))

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
                    mangled_rid,
                    θ;
                    RT,
                    ignore_reaction_ids,
                    ignore_metabolite_ids,
                )
            ) for (mangled_rid, (rid, rid_idx, mw)) in zip(reactions(smm)[col_idxs], kcat_idxs)
        ],
        n_reactions(smm),
    )

    #: inequality rhs
    Cp = coupling(smm.inner)
    M(θ) = [
        -1.0 * spdiagm(fill(1.0, num_reactions))
        1.0 * spdiagm(fill(1.0, num_reactions))
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

"""
$(TYPEDSIGNATURES)

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
