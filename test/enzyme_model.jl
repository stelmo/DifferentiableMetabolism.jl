using DifferentiableMetabolism, Test

import Symbolics as S
import AbstractFBCModels as AM
import ConstraintTrees as C
import COBREXA as X
import Tulip as T
import JuMP as J


Symbolics.@variables kcats_forward[1:4] kcats_backward[1:4] capacitylimitation

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
        upper_bound = 100.0,
        stoichiometry = Dict("m2" => 1.0),
    ),
    "r3" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m1" => -1.0, "m2" => -1.0, "m3" => 1.0),
    ),
    "r4" => AM.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m3" => -1.0, "m4" => 1.0),
    ),
    "r5" => AM.CanonicalModel.Reaction(
        lower_bound = -100.0,
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

reaction_isozymes = Dict(
    "r3" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[1], kcats_backward[1]),
    ),
    "r4" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[2], kcats_backward[2]),
        "iso2" =>
            ParameterIsozyme(Dict("g2" => 1), kcats_forward[3], kcats_backward[3]),
    ),
    "r5" => Dict(
        "iso1" => ParameterIsozyme(
            Dict("g3" => 1, "g4" => 2),
            kcats_forward[4],
            kcats_backward[4],
        ),
    ),
)

gene_molar_masses = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0, "g5" => 1.0)


# build differentiable model
m = COBREXA.fbc_model_constraints(model)
m += :enzymes^COBREXA.enzyme_variables(model)
m = COBREXA.add_enzyme_constraints!(m, reaction_isozymes)
m *=
    :total_proteome_bound^ConstraintTrees.Constraint(
        value = sum(
            m.enzymes[Symbol(gid)].value * gene_molar_masses[gid] for gid in AM.genes(model)
        ),
        bound = ParameterBetween(0.0, capacitylimitation),
    )

# substitute params into model
parameters = Dict(
    kcats_forward[1] => 1.0,
    kcats_forward[2] => 2.0,
    kcats_forward[3] => 3.0,
    kcats_forward[4] => 70.0,
    kcats_backward[1] => 1.0,
    kcats_backward[2] => 2.0,
    kcats_backward[3] => 3.0,
    kcats_backward[4] => 70.0,
    capacitylimitation => 0.5,
)

_x, _ν, _λ = optimized_constraints_with_parameters(
    m,
    parameters;
    objective = m.objective.value,
    optimizer = T.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
    duals = true,
)

@test isapprox(ec_solution.objective, 3.181818181753438, atol = 1e-3)
@test isapprox(ec_solution.enzymes.g4, 0.09090909090607537, atol = 1e-3)


kktfunc = kkt(m, m.objective.value)




