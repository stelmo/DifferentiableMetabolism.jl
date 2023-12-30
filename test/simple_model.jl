
import AbstractFBCModels as AM

# build simple model 
mets = Dict(
    "m1" => AM.CanonicalModel.Metabolite(),
    "m2" => AM.CanonicalModel.Metabolite(),
    "m3" => AM.CanonicalModel.Metabolite(),
    "m4" => AM.CanonicalModel.Metabolite(),
)

rxns = Dict(
    "r1" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m1" => 1.0),
    ),
    "r2" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 1.0,
        stoichiometry = Dict("m2" => 1.0),
    ),
    "r3" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m1" => -1.0, "m2" => -1.0, "m3" => 1.0),
    ),
    "r4" => AM.CanonicalModel.Reaction(
        lower_bound = -100.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m3" => -1.0, "m4" => 1.0),
    ),
    "r5" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m2" => -1.0, "m4" => 1.0),
    ),
    "r6" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m4" => -1.0),
        objective_coefficient = 1.0,
    ),
)

gs = Dict("g$i" => AM.CanonicalModel.Gene() for i = 1:5)

model = AM.CanonicalModel.Model(rxns, mets, gs)
