@testset "Differentiable GECKO LP" begin

    #: Set problem up
    model, protein_masses, protein_stoichiometry, reaction_kcats =
        import_core_model_and_data()

    model.reactions["EX_glc__D_e"].lb = -1000.0 # unconstraint because enzyme constraints take over
    total_protein_mass = 100 # mg/gdW

    remove_slow_isozymes!(model; reaction_kcats, protein_stoichiometry, protein_masses)

    obj_id = "BIOMASS_Ecoli_core_w_GAM§FOR"
    rxn_fluxes, prot_concens = gecko(
        model,
        CPLEX.Optimizer;
        objective_id = obj_id,
        protein_stoichiometry,
        protein_masses,
        reaction_kcats,
        total_protein_mass,
    )

    pruned_model = prune_model(model, rxn_fluxes; rtol = 1e-10)

    #: Differentiate an optimal solution

    kcat_rid_order = [
        rid for rid in reactions(pruned_model) if
        haskey(reaction_kcats, rid) && COBREXA._has_grr(pruned_model, rid)
    ]

    c, Ef, d, M, hf, reaction_map, protein_ids = differentiable_gecko_opt_problem(
        pruned_model;
        protein_stoichiometry,
        protein_masses,
        reaction_kcats,
        kcat_rid_order,
        ϵ = 1e-8,
    )

    θ = [
        [first(reaction_kcats[rid][1]) for rid in kcat_rid_order]
        total_protein_mass
    ]

    bid = reaction_map[obj_id]
    c[bid] = -1.0

    _, dx, _ = differentiate_LP(c, Ef, d, M, hf, θ, CPLEX.Optimizer)

    @test isapprox(sum(round.(dx, digits = 6)), 154.65145; atol = TEST_TOLERANCE)
end

@testset "Differentiable GECKO QP" begin
    using COBREXA, DifferentiableMetabolism, CPLEX, JSON, ForwardDiff

    #: Set problem up
    model, protein_masses, protein_stoichiometry, reaction_kcats =
        import_core_model_and_data()
    fluxdata = JSON.parsefile(joinpath("data", "fluxdata.json"))
    proteindata = JSON.parsefile(joinpath("data", "proteindata.json"))
    
    model.reactions["EX_glc__D_e"].lb = -1000.0 # unconstraint because enzyme constraints take over
    total_protein_mass = 100 # mg/gdW
    
    remove_slow_isozymes!(model; reaction_kcats, protein_stoichiometry, protein_masses)
    
    objective_id = "BIOMASS_Ecoli_core_w_GAM§FOR"
    rxn_fluxes, prot_concens = gecko(
        model,
        CPLEX.Optimizer;
        objective_id,
        protein_stoichiometry,
        protein_masses,
        reaction_kcats,
        total_protein_mass,
    )
    
    pruned_model = prune_model(model, rxn_fluxes; rtol = 1e-10)
    
    #: Differentiate an optimal solution
    
    kcat_rid_order = [
        rid for rid in reactions(pruned_model) if
        haskey(reaction_kcats, rid) && COBREXA._has_grr(pruned_model, rid)
    ]
    
    _, Ef, d, M, hf, reaction_map, protein_ids = differentiable_gecko_opt_problem(
        pruned_model;
        protein_stoichiometry,
        protein_masses,
        reaction_kcats,
        kcat_rid_order,
        ϵ = 1e-9,
    )
    
    θ = [
        [first(reaction_kcats[rid][1]) for rid in kcat_rid_order]
        total_protein_mass
    ]
    
    rids = [first(split(rid, "§")) for rid in COBREXA._order_id_to_idx_dict(reaction_map)]
    gids = protein_ids
    obs_v_dict = Dict(k => v[1] for (k, v) in fluxdata)
    obs_e_dict = Dict(k => v * 1e-9 for (k, v) in proteindata)
    Q, c, n = qp_objective_measured(
        rids,
        gids,
        obs_v_dict,
        obs_e_dict;
        vtol = 1e-3,
        etol = 1e-3,
        reg = 1e-1,
    )
    
    x, dx, obj = differentiate_QP(
        Q,
        c,
        n,
        Ef,
        d,
        M,
        hf,
        θ,
        CPLEX.Optimizer;
        modifications = [
            change_optimizer_attribute("CPXPARAM_Emphasis_Numerical", 1),
            COBREXA.silence
        ],
    )
    
    @test isapprox(sum(round.(dx, digits = 6)), -75.26632699999999; atol=TEST_TOLERANCE)
    @test isapprox(obj, 20.3157903402984; atol = TEST_TOLERANCE)
end
