function create_test_model()
    #=
    Implement the small model based on model found in the supplement in the
    original GECKO paper. This model is nice to troubleshoot with, because the
    stoich matrix is small. All the reactions are active by design and each reaction
    that has a kcat has only one grr.
    =#
    m = StandardModel("SmallModel")
    m1 = Metabolite("m1")
    m2 = Metabolite("m2")
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")
    m5 = Metabolite("m5")
    m6 = Metabolite("m6")

    @add_reactions! m begin
        "r1", nothing → m1, 0, 100
        "r2", nothing → m2, 0, 100
        "r3", m1 + m2 → m3, 0, 100
        "r4", m3 → m4 + m5, 0, 100
        "r5", m2 → m4 + m6, 0, 100
        "r6", m4 → nothing, 0, 100
        "biomass", m6 + m5 → nothing, 0, 100
    end

    gs = [Gene("g$i") for i = 1:4]

    m.reactions["biomass"].objective_coefficient = 1.0

    add_genes!(m, gs)
    add_metabolites!(m, [m1, m2, m3, m4, m5, m6])

    reaction_isozymes = Dict(
        "r3" => [Isozyme(Dict("g1" => 1), 10.0, 10.0)],
        "r4" => [Isozyme(Dict("g2" => 1, "g3" => 3), 30.0, 20.0)],
        "r5" => [Isozyme(Dict("g3" => 1, "g4" => 2), 70.0, 30.0)],
    )
    gene_product_bounds = Dict(
        "g1" => (0.0, 10.0),
        "g2" => (0.0, 10.0),
        "g3" => (0.0, 10.0),
        "g4" => (0.0, 10.0),
    )

    gene_product_molar_mass = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0)

    gene_product_mass_group_bound = Dict("uncategorized" => 1.0)

    return m,
    reaction_isozymes,
    gene_product_bounds,
    gene_product_molar_mass,
    gene_product_mass_group_bound
end
