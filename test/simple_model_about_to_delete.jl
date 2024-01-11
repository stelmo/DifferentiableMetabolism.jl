include("../../test/simple_model.jl") #hide

# ![simple_model](./assets/simple_model.svg)

model

# ## Add enzyme kinetic information

Symbolics.@variables kcats_forward[1:4] kcats_backward[1:4]

reaction_isozymes = Dict(
    "r3" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[1], kcats_backward[1]),
    ),
    "r4" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g2" => 1), kcats_forward[2], kcats_backward[2]),
        "iso2" =>
            ParameterIsozyme(Dict("g3" => 1), kcats_forward[3], kcats_backward[3]),
    ),
    "r5" => Dict(
        "iso1" => ParameterIsozyme(
            Dict("g4" => 1, "g5" => 2),
            kcats_forward[4],
            kcats_backward[4],
        ),
    ),
)

gene_molar_masses = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0, "g5" => 1.0)

# ## Add a parameterized capacity limitation

Symbolics.@variables capacitylimitation

# Build differentiable enzyme constrained model
m = COBREXA.fbc_model_constraints(model)
m += :enzymes^COBREXA.enzyme_variables(model)
m = COBREXA.add_enzyme_constraints!(m, reaction_isozymes)
m *=
    :total_proteome_bound^enzyme_capacity(
        m.enzymes,
        gene_molar_masses,
        AbstractFBCModels.genes(model),
        ParameterBetween(0, capacitylimitation),
    )
m.fluxes.r2.bound = ConstraintTrees.Between(0, 100) # remove typical FBA constraint on mass input

# substitute params into model
parameter_values = Dict(
    kcats_forward[1] => 1.0,
    kcats_forward[2] => 2.0,
    kcats_forward[3] => 3.0,
    kcats_forward[4] => 2.0,
    kcats_backward[1] => 1.0,
    kcats_backward[2] => 2.0,
    kcats_backward[3] => 3.0,
    kcats_backward[4] => 2.0,
    capacitylimitation => 0.5,
)

ec_solution, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    m,
    parameter_values;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

# @test isapprox(ec_solution.objective, 3.1818181867946134, atol = TEST_TOLERANCE)
# @test isapprox(ec_solution.enzymes.g4, 0.09090909089739453, atol = TEST_TOLERANCE)
# @test isapprox(ec_solution.:total_proteome_bound, 0.5, atol = TEST_TOLERANCE)

ec_solution.fluxes
ec_solution.enzymes

sens = differentiate(
    m,
    m.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    [capacitylimitation; kcats_forward; kcats_backward],
)

# make solution unique

reaction_isozymes = Dict(
    "r3" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[1], kcats_backward[1]),
    ),
    "r4" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g2" => 1), kcats_forward[2], kcats_backward[2]),
        "iso2" =>
            ParameterIsozyme(Dict("g3" => 1), kcats_forward[3], kcats_backward[3]),
    ),
    "r5" => Dict(
        "iso1" => ParameterIsozyme(
            Dict("g4" => 1, "g5" => 2),
            kcats_forward[4],
            kcats_backward[4],
        ),
    ),
)

delete!(reaction_isozymes["r4"], "iso2")

m = COBREXA.fbc_model_constraints(model)
m += :enzymes^COBREXA.enzyme_variables(model)
m = COBREXA.add_enzyme_constraints!(m, reaction_isozymes)
m *=
    :total_proteome_bound^enzyme_capacity(
        m.enzymes,
        gene_molar_masses,
        AbstractFBCModels.genes(model),
        ParameterBetween(0, capacitylimitation),
    )
m.fluxes.r2.bound = ConstraintTrees.Between(0, 100) # remove typical FBA constraint on mass input


ec_solution, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    m,
    parameter_values;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution.fluxes
ec_solution.enzymes

sens = differentiate(
    m,
    m.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    [capacitylimitation; kcats_forward; kcats_backward],
)


include("../../test/simple_model.jl") #hide

# ![simple_model](./assets/simple_model.svg)

model

# ## Add enzyme kinetic information

Symbolics.@variables kcats_forward[1:4] kcats_backward[1:4]

reaction_isozymes = Dict(
    "r3" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g1" => 1), kcats_forward[1], kcats_backward[1]),
    ),
    "r4" => Dict(
        "iso1" =>
            ParameterIsozyme(Dict("g2" => 1), kcats_forward[2], kcats_backward[2]),
        "iso2" =>
            ParameterIsozyme(Dict("g3" => 1), kcats_forward[3], kcats_backward[3]),
    ),
    "r5" => Dict(
        "iso1" => ParameterIsozyme(
            Dict("g4" => 1, "g5" => 2),
            kcats_forward[4],
            kcats_backward[4],
        ),
    ),
)

gene_molar_masses = Dict("g1" => 1.0, "g2" => 2.0, "g3" => 3.0, "g4" => 4.0, "g5" => 1.0)

# ## Add a parameterized capacity limitation

Symbolics.@variables capacitylimitation

# Build differentiable enzyme constrained model
m = COBREXA.fbc_model_constraints(model)
m += :enzymes^COBREXA.enzyme_variables(model)
m = COBREXA.add_enzyme_constraints!(m, reaction_isozymes)
m *=
    :total_proteome_bound^enzyme_capacity(
        m.enzymes,
        gene_molar_masses,
        AbstractFBCModels.genes(model),
        ParameterBetween(0, capacitylimitation),
    )
m.fluxes.r2.bound = ConstraintTrees.Between(0, 100) # remove typical FBA constraint on mass input

# substitute params into model
parameter_values = Dict(
    kcats_forward[1] => 1.0,
    kcats_forward[2] => 2.0,
    kcats_forward[3] => 3.0,
    kcats_forward[4] => 2.0,
    kcats_backward[1] => 1.0,
    kcats_backward[2] => 2.0,
    kcats_backward[3] => 3.0,
    kcats_backward[4] => 2.0,
    capacitylimitation => 0.5,
)

ec_solution, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    m,
    parameter_values;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution.objective

# add small quadratic weight
m.objective = ConstraintTrees.Constraint(
    value = 0.001 * (
        sum(x.value * x.value for x in values(m.fluxes)) +
        sum(x.value * x.value for x in values(m.enzymes))
    ) - m.objective.value,
    bound = nothing,
)

using Gurobi
ec_solution, x_vals, eq_dual_vals, ineq_dual_vals = optimized_constraints_with_parameters(
    m,
    parameter_values;
    objective = m.objective.value,
    optimizer = Gurobi.Optimizer,
    sense = COBREXA.Minimal,
)

# @test isapprox(ec_solution.objective, 3.1818181867946134, atol = TEST_TOLERANCE)
# @test isapprox(ec_solution.enzymes.g4, 0.09090909089739453, atol = TEST_TOLERANCE)
# @test isapprox(ec_solution.:total_proteome_bound, 0.5, atol = TEST_TOLERANCE)

ec_solution.fluxes
ec_solution.enzymes

sens = differentiate(
    m,
    m.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    [capacitylimitation; kcats_forward; kcats_backward],
)

objective = m.objective.value
parameters = [capacitylimitation; kcats_forward; kcats_backward]
