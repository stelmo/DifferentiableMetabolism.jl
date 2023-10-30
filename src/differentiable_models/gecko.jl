"""
$(TYPEDSIGNATURES)

Construct a [`DifferentiableModel`](@ref) from a [`COBREXA.GeckoModel`](@ref).
Each variable in `gm` is differentiated with respect to the kcats in the
dictionary `rid_enzyme`, which is a dictionary mapping reaction ids to
[`Enzyme`](@ref)s. Enzyme constraints are only taken with respect to the entries
of `rid_enzyme`. Optionally, incorporate thermodynamic and saturation constraints through
`rid_dg0` and `rid_km`. If either thermodynamic or saturation (or both) constraints 
are added, then metabolite concentrations, `mid_concentration`, need to be 
supplied as well.

Note, thermodynamic parameters require special attention. To ignore some
reactions when calculating the thermodynamic factor, include them in
`ignore_reaction_ids`. The units used for the Gibbs free energy are kJ/mol, and
concentrations should be in molar. Importantly, the standard Gibbs free energy
of reaction is assumed to be given *in the forward direction of the reaction* in
the model.

Internally, `ϵ`, `atol`, and `digits`, are forwarded to
[`_remove_lin_dep_rows`](@ref). Optionally, scale the equality constraints with
`scale_equality`, see [`scaling_factor`](@ref) for more information.
"""
function with_parameters(
    gm::GeckoModel,
    rid_enzyme::Dict{String,Enzyme};
    rid_dg0 = Dict{String,Float64}(),
    rid_km = Dict{String,Dict{String,Float64}}(),
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
    if isempty(rid_dg0) && isempty(rid_km)
        # no concentration parameters
        param_ids = "k#" .* collect(keys(rid_enzyme))
        θ = [x.kcat for x in values(rid_enzyme)]
    else
        # has concentration parameters
        isempty(mid_concentration) && throw(error("Missing metabolite concentrations."))
        
        param_ids = [
            "k#" .* collect(keys(rid_enzyme))
            "c#" .* metabolites(gm.inner)
        ]
        θ = [
            [x.kcat for x in values(rid_enzyme)]
            [mid_concentration[mid] for mid in metabolites(gm.inner)]
        ]    
    end

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
    ignore_metabolite_ids = [],
)
    #: get irreverible stoichiometric matrix from model
    irrev_S = stoichiometry(gm.inner) * COBREXA._gecko_reaction_column_reactions(gm)

    #: make full rank
    S = DifferentiableMetabolism._remove_lin_dep_rows(irrev_S; ϵ, digits, atol)

    #: size of resultant model
    num_reactions = n_reactions(gm) - n_genes(gm)
    num_genes = n_genes(gm)
    num_metabolites = size(S, 1) # take into account removed linear dependencies
    num_vars = n_reactions(gm)

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
                    mangled_rid,
                    θ;
                    RT,
                    ignore_reaction_ids,
                    ignore_metabolite_ids,
                ) * 
                DifferentiableMetabolism._saturation(
                    gm,
                    rid_enzyme,
                    rid_km,
                    rid,
                    mangled_rid,
                    θ;
                    ignore_reaction_ids,
                    ignore_metabolite_ids,
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

    return c, E, d, M, h, reactions(gm)
end

"""
$(TYPEDSIGNATURES)

Helper function to add a column into the enzyme stoichiometric matrix
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
$(TYPEDSIGNATURES)

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