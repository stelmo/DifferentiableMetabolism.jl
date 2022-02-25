function import_core_model_and_data()
    # Import test models 
    model = load_model(StandardModel, joinpath("data", "e_coli_core_modified_model.json"))
    protein_masses = JSON.parsefile(joinpath("data", "protein_masses.json"))
    protein_stoichiometry = JSON.parsefile(joinpath("data", "protein_stoichiometry.json"))
    reaction_kcats = JSON.parsefile(joinpath("data", "reaction_kcats.json"))

    return model, protein_masses, protein_stoichiometry, reaction_kcats
end
