@testset "Differentiable SMOMENT" begin
    #: Set problem up
    model, protein_masses, reaction_protein_stoichiometry, reaction_kcats =
        import_core_model_and_data()

    model.reactions["EX_glc__D_e"].lb = -1000.0 # unconstraint because enzyme constraints take over
    total_protein_mass = 100 # mg/gdW

    remove_slow_isozymes!(
        model;
        reaction_kcats,
        reaction_protein_stoichiometry,
        protein_masses,
    )

    smm = SMomentModel(
        model;
        reaction_kcats,
        reaction_protein_stoichiometry,
        protein_masses,
        total_protein_mass, # mg/gdW
    )

    rxn_fluxes = flux_balance_analysis_dict(smm, CPLEX.Optimizer)

    smm.smodel = prune_model(model, rxn_fluxes; rtol = 1e-10)

    #: Differentiate an optimal solution
    optimizer = CPLEX.Optimizer
    res = differentiate_smoment(smm, optimizer)

    @test isapprox(
        sum(round.(res.dx, digits = 6)),
        95.23119400000002;
        atol = TEST_TOLERANCE,
    )
end
