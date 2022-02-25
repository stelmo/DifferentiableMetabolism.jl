@testset "Differentiable Thermodynamic SMOMENT" begin

    #: Set problem up
    model, protein_masses, protein_stoichiometry, reaction_kcats =
        import_core_model_and_data()
    
    #: Load ΔG⁰s for all reactions in model that have thermodata
    reaction_dg0s = Dict(
        k => float(v) for
        (k, v) in JSON.parsefile(joinpath("data", "e_coli_core_dgs.json"))
    )
    
    #: find a thermodynamically consistent solution
    loopless_sol = flux_balance_analysis_dict(
        model,
        CPLEX.Optimizer;
        modifications = [add_loopless_constraints()],
    )
    
    #: Find concentration profile that explains the loopless solution
    mmdf = max_min_driving_force(
        model,
        reaction_dg0s,
        CPLEX.Optimizer;
        flux_solution = loopless_sol,
        proton_ids = ["h_c", "h_e"],
        water_ids = ["h2o_c", "h2o_e"],
        concentration_ratios = Dict(
            ("atp_c", "adp_c") => 10.0,
            ("nadh_c", "nad_c") => 0.13,
            ("nadph_c", "nadp_c") => 1.3,
        ),
        concentration_lb = 1e-9,
        concentration_ub = 100e-3,
        modifications = [],
        ignore_reaction_ids = ["H2Ot"],
    )
    metabolite_concentrations = mmdf.concentrations
    
    #: change model to incorporate flux directions (via pruning)
    remove_slow_isozymes!(model; reaction_kcats, protein_stoichiometry, protein_masses)
    model = prune_model(model, loopless_sol; rtol = 1e-10)
    
    #: Thermokinetic SMOMENT
    total_protein_mass = 100.0 # mg/gDW
    model.reactions["EX_glc__D_e"].lb = -1000.0 #! unconstrain otherwise bound will be hit
    objective_id = "BIOMASS_Ecoli_core_w_GAM§FOR"
    
    rxn_fluxes = thermodynamic_smoment(
        model,
        CPLEX.Optimizer;
        objective_id,
        metabolite_concentrations,
        protein_stoichiometry,
        protein_masses,
        reaction_kcats,
        reaction_dg0s,
        total_protein_mass,
        ignore_reaction_ids = ["H2Ot"],
    )
    
    #: prune model
    model = prune_model(model, rxn_fluxes; rtol = 1e-10)
    
    #: Differentiate pruned model
    kcat_rid_order = [
        rid for rid in reactions(model) if
        haskey(reaction_kcats, rid) && COBREXA._has_grr(model, rid)
    ]
    
    met_conc_order =
        [mid for mid in metabolites(model) if haskey(metabolite_concentrations, mid)] 
    #! assumption is that every metabolite here is involved in some reaction
    
    c, Ef, d, M, hf, reaction_map, metabolite_map = differentiable_thermokinetic_smoment(
        model;
        protein_stoichiometry,
        protein_masses,
        kcat_rid_order,
        met_conc_order,
        reaction_dg0s,
        reaction_kcats,
        ϵ = 1e-8,
        ignore_reaction_ids = ["H2Ot"],
    )
    
    θ = [
        [first(reaction_kcats[rid][1]) for rid in kcat_rid_order]
        [metabolite_concentrations[mid] for mid in met_conc_order]
        total_protein_mass
    ]
    
    bid = reaction_map[objective_id]
    c[bid] = -1.0
    
    x, dx, obj = differentiate_LP(c, Ef, d, M, hf, θ, CPLEX.Optimizer)
    
    @test isapprox(sum(round.(dx, digits = 6)), 105.96230700000001; atol = TEST_TOLERANCE)
    @test isapprox(obj, -0.5329382139620737; atol = TEST_TOLERANCE)
end
