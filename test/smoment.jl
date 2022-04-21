@testset "Differentiable SMOMENT" begin
    #: Set problem up 
    model, rid_isozymes = import_core_model_and_data();
    change_bound!(model, "EX_glc__D_e"; lower = -1000.0)
    total_protein_mass = 100 # mg/gdW
    
    remove_slow_isozymes!(
        model,
        rid_isozymes,
    )
    smm = SMomentModel(
        model;
        rid_isozymes,
        enzyme_capacity = total_protein_mass,
    )

    rxn_fluxes = flux_balance_analysis_dict(smm, CPLEX.Optimizer)

    pruned_model = prune_model(model, rxn_fluxes)
    
    #: remove isozymes no longer in model 
    delete!.(
        Ref(rid_isozymes), 
        filter(x->!haskey(pruned_model.reactions, x), keys(rid_isozymes)),
    );
    
    simplified_smm = SMomentModel(
        pruned_model;
        rid_isozymes,
        enzyme_capacity = total_protein_mass,
    ) 
    rxn_fluxes = flux_balance_analysis_dict(simplified_smm, CPLEX.Optimizer)

    rid_enzyme = isozyme_to_enzyme(model, rid_isozymes) 

    #: differentiate model 
    diffmodel = with_parameters(
        simplified_smm, 
        rid_enzyme,
    ) 

    x, dx = differentiate(
        diffmodel,
        CPLEX.Optimizer,
    )

    @test isapprox(
        sum(round.(dx, digits = 6)),
        95.23119400000002;
        atol = TEST_TOLERANCE,
    )
end
