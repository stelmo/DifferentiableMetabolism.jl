using DifferentiableMetabolism

import Symbolics as S
import AbstractFBCModels as A
import ConstraintTrees as C
import COBREXA as X
import Tulip as T

S.@variables kcats_forward[1:4] kcats_backward[1:4] capacitylimitation

# build simple model 
mets = Dict(
    "m1" => A.CanonicalModel.Metabolite(),
    "m2" => A.CanonicalModel.Metabolite(),
    "m3" => A.CanonicalModel.Metabolite(),
    "m4" => A.CanonicalModel.Metabolite(),
)

rxns = Dict(
    "r1" => A.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m1" => 1.0),
    ),
    "r2" => A.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m2" => 1.0),
    ),
    "r3" => A.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m1" => -1.0, "m2" => -1.0, "m3" => 1.0),
    ),
    "r4" => A.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m3" => -1.0, "m4" => 1.0),
    ),
    "r5" => A.CanonicalModel.Reaction(
        lower_bound = -100.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m2" => -1.0, "m4" => 1.0),
    ),
    "r6" => A.CanonicalModel.Reaction(
        lower_bound = 0.0,
        upper_bound = 100.0,
        stoichiometry = Dict("m4" => -1.0),
        objective_coefficient = 1.0,
    ),
)

gs = Dict("g$i" => A.CanonicalModel.Gene() for i = 1:5)

model = A.CanonicalModel.Model(rxns, mets, gs)

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
m = X.fbc_model_constraints(model)
m += :enzymes^X.enzyme_variables(model)
m = X.add_enzyme_constraints!(m, reaction_isozymes)
m *= :total_proteome_bound^C.Constraint(
    value = sum(
        m.enzymes[Symbol(gid)].value * gene_molar_masses[gid] for gid in A.genes(model)
    ),
    bound = ParameterBetween(0.0, capacitylimitation),
)

# full symbolic model
S.@variables x[1:C.var_count(m)]
sm = C.constraint_values(m, collect(x))

params = Dict(
    kcats_forward[1] => 1.0,
    kcats_forward[2] => 2.0,
    kcats_forward[3] => 3.0,
    kcats_forward[4] => 70.0,
    kcats_backward[1] => 1.0,
    kcats_backward[2] => 2.0,
    kcats_backward[3] => 3.0,
    kcats_backward[4]  => 70.0,
    capacitylimitation => 0.5,
)
pm = S.substitute(m, params) #  type of tree is Any which does not work downstream

# solve the model
ec_solution = X.optimized_constraints(
    pm;
    objective = pm.objective.value,
    optimizer = T.Optimizer,
    modifications = [X.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)




S.@variables p[1:4]
m = :variables^C.variables(keys = [:x, :y])
m *= :param_constraints^C.ConstraintTree(
    :c1 => C.Constraint(
        value = p[1] * m.variables.x.value + m.variables.y.value, 
        bound = ParameterBetween(0.0, p[2]),
    ),
    :c2 => C.Constraint(
        value = m.variables.x.value + p[3] * m.variables.y.value, 
        bound = ParameterEqualTo(p[4]),
    ),
)

params = Dict(
    p[1] => 1.0,
    p[2] => 2.0,
    p[3] => 3.0,
    p[4] => 4.0,
)
pm = C.substitute(m, params) #  type of tree is Any which does not work downstream
S.substitute(m.param_constraints.c1, params) # this works
