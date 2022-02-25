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

    c, Ef, d, M, hf, reaction_map, protein_ids = differentiable_gecko(
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

end
