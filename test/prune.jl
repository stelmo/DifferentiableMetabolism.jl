@testset "Prune model" begin

    m = StandardModel("SmallModel")
    m1 = Metabolite("m1")
    m2 = Metabolite("m2")
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")
    m5 = Metabolite("m5")

    @add_reactions! m begin
        "r1", nothing → m1, 0, 100
        "r2", nothing → m2, 0, 100
        "r3", m1 + m2 → m3, 0, 100
        "r4", m3 → m4, 0, 100
        "r5", m3 → m5, 0, 100
        "r6", m5 → m4, 0, 100
        "r7", m4 → nothing, 0, 100
    end

    gs = [Gene("g$i") for i = 1:4]

    m.reactions["r3"].grr = [["g1"], ["g2"]]
    m.reactions["r4"].grr = [["g3", "g4"]]
    m.reactions["r7"].objective_coefficient = 1.0

    add_genes!(m, gs)
    add_metabolites!(m, [m1, m2, m3, m4, m5])

    reaction_fluxes = Dict(
        "r1" => 1.0,
        "r2" => 1.0,
        "r3" => -1.0,
        "r4" => 1.0,
        "r5" => -0.01,
        "r6" => 0.01,
        "r7" => 1.0,
    )

    pruned_model = prune_model(m, reaction_fluxes; atol = 1e-2)

    @test !haskey(pruned_model.reactions, "r5")
    @test !haskey(pruned_model.reactions, "r6")
    @test haskey(pruned_model.reactions, "r1")
    @test !isempty(genes(pruned_model))
    @test "m1" in metabolites(pruned_model)
    @test "m5" ∉ metabolites(pruned_model)
end
