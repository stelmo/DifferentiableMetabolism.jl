@testset "Differentiable Thermodynamic SMOMENT" begin

    #: Set problem up
    model, protein_masses, reaction_protein_stoichiometry, reaction_kcats =
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
    remove_slow_isozymes!(
        model;
        reaction_kcats,
        reaction_protein_stoichiometry,
        protein_masses,
    )

    model = prune_model(model, loopless_sol; rtol = 1e-10)
    total_protein_mass = 100.0 # mg/gDW
    model.reactions["EX_glc__D_e"].lb = -1000.0 #! unconstrain otherwise bound will be hit

    smm = SMomentModel(
        model;
        reaction_kcats,
        reaction_protein_stoichiometry,
        protein_masses,
        total_protein_mass, # mg/gdW
    )

    #: Differentiate pruned model
    res = differentiate_thermodynamic_smoment(
        smm,
        CPLEX.Optimizer;
        reaction_dg0s,
        metabolite_concentrations,
        ignore_reaction_ids = ["H2Ot"],
        scale_input = false,
        scale_output = true,
    )
    rxn_fluxes = COBREXA._map_irrev_to_rev_ids(smm.smomentdata.reaction_map, res.x)

    @test isapprox(
        sum(round.(res.dx, digits = 6)),
        105.96230700000001;
        atol = TEST_TOLERANCE,
    )
    @test isapprox(
        rxn_fluxes["BIOMASS_Ecoli_core_w_GAM"],
        0.5329382139620737;
        atol = TEST_TOLERANCE,
    )
end
