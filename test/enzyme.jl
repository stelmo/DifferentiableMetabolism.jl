@testset "Enzyme type and conversion" begin
    # create dummy isozyme
    isozyme = Isozyme(Dict("g1" => 1, "g2" => 3), 10.0, 2.0)

    # create dummy molar masses
    gene_product_mass_lookup = Dict("g1" => 5.0, "g2" => 7.0)

    # convert isozyme to enzyme with forward kcat
    enz_for = isozyme_to_enzyme(isozyme, gene_product_mass_lookup; direction = :forward)

    @test enz_for.kcat == 10.0
    @test enz_for.gene_product_count["g1"] == 1
    @test enz_for.gene_product_count["g2"] == 3
    @test enz_for.gene_product_mass["g1"] == 5.0
    @test enz_for.gene_product_mass["g2"] == 7.0

    # convert isozyme to enzyme with reverse kcat
    enz_rev = isozyme_to_enzyme(isozyme, gene_product_mass_lookup; direction = :reverse)

    @test enz_rev.kcat == 2.0
end
