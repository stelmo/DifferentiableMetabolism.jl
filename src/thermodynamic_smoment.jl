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

    c, E, d, M, h, var_ids = differentiable_thermodynamic_smoment_opt_problem(
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
function differentiable_thermodynamic_smoment_opt_problem(
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
    kcat_original_rid_order = String[]
    met_order = Vector{Vector{String}}()
    nus_vec = Vector{Vector{Float64}}()
    col_idxs = Int[]
    mws = Float64[]
    for (col_idx, rid) in enumerate(reactions(smm))
        original_rid = string(first(split(rid, "#")))

        # skip these entries
        !haskey(rid_enzyme, original_rid) && continue
        original_rid in ignore_reaction_ids && continue

        # kinetic info
        mw = sum([
            pstoich * rid_enzyme[original_rid].gene_product_mass[gid] for
            (gid, pstoich) in rid_enzyme[original_rid].gene_product_count
        ])
        push!(kcat_original_rid_order, original_rid)
        push!(col_idxs, col_idx)
        push!(mws, mw)

        # thermo info
        rs = reaction_stoichiometry(smm.inner, original_rid)
        if haskey(rid_dg0, original_rid)
            push!(met_order, collect(keys(rs)))
            push!(nus_vec, collect(values(rs)))
        else
            push!(met_order, collect(keys(rs)))
            push!(nus_vec, Float64[]) # acts like a flag
        end
    end

    kcat_idxs = [
        (rid, first(indexin([rid], collect(keys(rid_enzyme)))), mw) for
        (rid, mw) in zip(kcat_original_rid_order, mws)
    ]

    mid_idxs = [
        (nus, Int.(indexin(mids, metabolites(smm))) .+ length(rid_enzyme)) for
        (nus, mids) in zip(nus_vec, met_order)
    ]

    kcat_thermo_coupling(θ) = sparsevec(
        col_idxs,
        [
            mw / (
                (θ[rid_idx]) *
                DifferentiableMetabolism._dg(rid, nus, θ[mconcs_idxs], rid_dg0, RT)
            ) for ((rid, rid_idx, mw), (nus, mconcs_idxs)) in zip(kcat_idxs, mid_idxs)
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

"""
    $(TYPEDSIGNATURES)

Helper function to assign the thermodynamic driving to a reaction.
"""
function _dg(rid, nus, mconcs, rid_dg0s, RT)
    # return 1.0
    isempty(nus) && return 1.0

    dg_val = rid_dg0s[rid] + RT * sum(nu * log(mconc) for (nu, mconc) in zip(nus, mconcs))
    return 1.0 - exp(dg_val / RT)
end
