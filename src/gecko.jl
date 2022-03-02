"""


Notes:
1) `model` must be pruned, i.e. only active reactions, metabolites and genes.
Singularity issues will arise if not.
2) Each reaction has only one active isozyme and one kcat assigned to that
isozyme.
3) Rescale some columns to prevent numerical issues from arising
"""
function differentiable_gecko_opt_problem(
    model::StandardModel;
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    lb_protein_measurements = Dict(),
    ub_protein_measurements = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    kcat_rid_order = [],
    ϵ = 1e-8,
    stoich_digits_round = 8,
)
    S_rank_issues, lb_fluxes, ub_fluxes, reaction_map, _ =
        COBREXA._build_irreversible_stoichiometric_matrix(model)

    #: make full rank
    S = round.(_remove_lin_dep_rows(S_rank_issues; ϵ), digits = stoich_digits_round)

    #: find all gene products that have kcats associated with them
    protein_ids = COBREXA._get_proteins_with_kcats(model, reaction_kcats)

    #: size of resultant model
    n_reactions = size(S, 2)
    n_proteins = length(protein_ids)
    n_metabolites = size(S, 1)
    n_vars = n_reactions + n_proteins

    #: equality lhs
    E_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeff_tuple = Vector{Tuple{Float64,String}}(),
    )

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        if contains(rid, "§ARM")
            @warn "Isozyme found, unexpected!"
            continue
        end
        !haskey(reaction_kcats, original_rid) && continue

        # add all entries to column of matrix
        #! no isozymes by assumption and all unidirectional
        _add_enzyme_variable_as_function(
            model,
            original_rid,
            protein_stoichiometry,
            E_components,
            col_idx,
            protein_ids,
        )
    end

    kcat_pst_idxs = [
        (pst, first(indexin([rid], kcat_rid_order))) for
        (pst, rid) in E_components.coeff_tuple
    ]
    Se(θ) = sparse(
        E_components.row_idxs,
        E_components.col_idxs,
        [pst / θ[x] for (pst, x) in kcat_pst_idxs],
        n_proteins,
        n_reactions,
    )

    Ef(θ) = Array([
        S zeros(n_metabolites, n_proteins)
        Se(θ) I(n_proteins)
    ])

    #: equality rhs
    d = zeros(n_metabolites + n_proteins) #TODO handle fixed variables

    #: need to set objective reaction outside
    c = zeros(n_vars)

    #: inequality constraints
    M, hf = _gecko_build_inequality_constraints_as_functions(
        lb_protein_measurements,
        ub_protein_measurements,
        protein_ids,
        protein_masses,
        n_reactions,
        n_proteins,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

    opt = (c = c, Ef = Ef, d = d, M = M, hf = hf)
    opt_struct = (
        reaction_map = reaction_map,
        protein_ids = protein_ids,
        kcat_rid_order = kcat_rid_order,
        row_idxs = E_components.row_idxs,
        col_idxs = E_components.col_idxs,
        vals = kcat_pst_idxs,
        n_proteins = n_proteins,
        n_reactions = n_reactions,
        n_metabolites = n_metabolites,
    )

    return opt, opt_struct
end

"""
    _add_enzyme_variable_as_function(model, rid, original_rid, protein_stoichiometry, reaction_kcats, E_components, col_idx, protein_ids)

Helper function to add an column into the enzyme stoichiometric matrix parametrically. Note, 
assumes only one grr.
"""
function _add_enzyme_variable_as_function(
    model,
    original_rid,
    protein_stoichiometry,
    E_components,
    col_idx,
    protein_ids,
)
    grr = reaction_gene_association(model, original_rid)[1] #! assumes both kcats are the same (appropriate one duplicated)
    pstoich = protein_stoichiometry[original_rid][1] #! requires pruned protein_stoichiometry
    for (idx, pid) in enumerate(grr)
        push!(E_components.row_idxs, first(indexin([pid], protein_ids)))
        push!(E_components.col_idxs, col_idx)
        push!(E_components.coeff_tuple, (-pstoich[idx], original_rid))
    end
end

"""
    _gecko_build_inequality_constraints_as_functions(
        lb_protein_measurements,
        ub_protein_measurements,
        protein_ids,
        protein_masses,
        n_reactions,
        n_proteins,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

Helper function to build inequality constraints. Returns M and h, but h is a
function taking a vector, θ (M is a matrix). The last element in θ is taken as
the total protein content and is the only part of θ that is used.
"""
function _gecko_build_inequality_constraints_as_functions(
    lb_protein_measurements,
    ub_protein_measurements,
    protein_ids,
    protein_masses,
    n_reactions,
    n_proteins,
    lb_flux_measurements,
    ub_flux_measurements,
    lb_fluxes,
    ub_fluxes,
    reaction_map,
)
    #: inequality lhs
    mw_proteins = [protein_masses[pid] for pid in protein_ids]
    M = Array(
        [
            -I(n_reactions) zeros(n_reactions, n_proteins)
            I(n_reactions) zeros(n_reactions, n_proteins)
            zeros(n_proteins, n_reactions) -I(n_proteins)
            zeros(n_proteins, n_reactions) I(n_proteins)
            zeros(1, n_reactions) mw_proteins'
        ],
    )

    #: inequality rhs
    for original_rid in keys(lb_flux_measurements) # only constrain if measurement available
        lb = lb_flux_measurements[original_rid]
        ub = ub_flux_measurements[original_rid]
        rids = [rid for rid in keys(reaction_map) if startswith(rid, original_rid)]

        any(contains(x, "§ISO") for x in rids) &&
            @warn("Isozyme detected, unexpected behaviour!")

        if lb > 0 # forward only
            for rid in rids
                contains(rid, "§REV") && (ub_fluxes[reaction_map[rid]] = 0.0)
                contains(rid, "§FOR") &&
                    (ub_fluxes[reaction_map[rid]] = ub; lb_fluxes[reaction_map[rid]] = lb)
            end
        elseif ub < 0 # reverse only
            for rid in rids
                contains(rid, "§FOR") && (ub_fluxes[reaction_map[rid]] = 0.0)
                contains(rid, "§REV") &&
                    (ub_fluxes[reaction_map[rid]] = -lb; lb_fluxes[reaction_map[rid]] = -ub)
            end
        else # measurement does not rule our reversibility
            for rid in rids
                contains(rid, "§FOR") &&
                    (ub_fluxes[reaction_map[rid]] = ub; lb_fluxes[reaction_map[rid]] = 0)
                contains(rid, "§REV") &&
                    (ub_fluxes[reaction_map[rid]] = -lb; lb_fluxes[reaction_map[rid]] = 0)
            end
        end
    end

    lb_proteins = [
        haskey(lb_protein_measurements, pid) ? lb_protein_measurements[pid] : 0.0 for
        pid in protein_ids
    ]
    ub_proteins = [
        haskey(ub_protein_measurements, pid) ? ub_protein_measurements[pid] : 1000.0 for
        pid in protein_ids
    ]

    hf(θ) = Array([-lb_fluxes; ub_fluxes; -lb_proteins; ub_proteins; θ[end]])

    return M, hf
end
