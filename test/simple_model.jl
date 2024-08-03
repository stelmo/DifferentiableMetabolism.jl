
import AbstractFBCModels
const AC = AbstractFBCModels.CanonicalModel

rs = [
    AC.Reaction(;
        name = "r1",
        lower_bound = -2,
        upper_bound = 0,
        stoichiometry = Dict("m1" => -1.0),
    )
    AC.Reaction(;
        name = "r2",
        lower_bound = -1,
        upper_bound = 0,
        stoichiometry = Dict("m2" => -1),
    )
    AC.Reaction(;
        name = "r3",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m1" => -1, "m2" => -1, "m3" => 1),
        gene_association_dnf = [["g1"]],
    )
    AC.Reaction(;
        name = "r4",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m3" => -1, "m4" => 1),
        gene_association_dnf = [["g2"], ["g3"]],
    )
    AC.Reaction(;
        name = "r5",
        lower_bound = -1000,
        upper_bound = 1000,
        stoichiometry = Dict("m2" => -1, "m4" => 1),
        gene_association_dnf = [["g4", "g5"]],
    )
    AC.Reaction(;
        name = "r6",
        lower_bound = 0,
        upper_bound = 1000,
        stoichiometry = Dict("m4" => -1),
        objective_coefficient = 1.0,
    )
]

ms = [
    AC.Metabolite(; name = "m1")
    AC.Metabolite(; name = "m2")
    AC.Metabolite(; name = "m3")
    AC.Metabolite(; name = "m4")
]

gs = [
    AC.Gene(; name = "g1")
    AC.Gene(; name = "g2")
    AC.Gene(; name = "g3")
    AC.Gene(; name = "g4")
    AC.Gene(; name = "g5")
]

model = AC.Model(
    reactions = Dict(x.name => x for x in rs),
    metabolites = Dict(x.name => x for x in ms),
    genes = Dict(x.name => x for x in gs),
)
