
# # Parametric constraint-based metabolic models

using DifferentiableMetabolism

using Symbolics
using ConstraintTrees
using COBREXA
using Tulip

# ## Load and solve a simple model

# The code used to construct the model is located in `test/simple_model.jl`, but
# it is not shown here for brevity. Below is a visualization of the model.

include(joinpath("..", "test", "simple_model.jl")) #hide

# ![simple_model](./assets/simple_model.svg)

model

# Build a basic ConstraintTree model without parameters
m = COBREXA.fbc_model_constraints(model)

# Solve normally
base_model = COBREXA.optimized_constraints(
    m;
    optimizer = Tulip.Optimizer,
    objective = m.objective.value,
)
base_model.fluxes

@test isapprox(base_model.objective, 1.0; atol = TEST_TOLERANCE)

# ## Add parameters to the model

# Make bound of r2 and mass balance of m3 parameters
Symbolics.@variables r2bound m3bound

m.fluxes.r2 = ConstraintTrees.Constraint(m.fluxes.r2.value, ParameterBetween(0, r2bound))

m.flux_stoichiometry.m3 =
    ConstraintTrees.Constraint(m.flux_stoichiometry.m3.value, ParameterEqualTo(m3bound))

# # add parametric constraints
Symbolics.@variables p[1:4]

m *=
    :linparam^ConstraintTrees.Constraint(
        value = p[1] * m.fluxes.r1.value + p[2] * m.fluxes.r2.value,
        bound = ParameterBetween(0, p[3]),
    )

# substitute params into model
parameter_substitutions = Dict(
    r2bound => 4.0,
    m3bound => 0.1, # lose some mass here
    p[1] => 1.0,
    p[2] => 1.0,
    p[3] => 4.0,
    p[4] => 1.0,
)

m_noparams = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
)
m_noparams.fluxes

@test isapprox(m_noparams.objective, 3.9; atol = TEST_TOLERANCE)

# ## Change the parameters and re-solve

# substitute params into model
parameter_substitutions[m3bound] = 0.0


m_noparams = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
)
m_noparams.fluxes

@test isapprox(m_noparams.objective, 4.0; atol = TEST_TOLERANCE)
