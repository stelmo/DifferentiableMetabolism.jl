```@meta
EditURL = "1-parametric-models.jl"
```

# Parametric constraint-based metabolic models

````@example 1-parametric-models
using DifferentiableMetabolism

using Symbolics
using ConstraintTrees
using COBREXA
using Tulip
using Clarabel
````

## Load and solve a simple model

The code used to construct the model is located in `test/simple_model.jl`, but
it is not shown here for brevity. Below is a visualization of the model.

````@example 1-parametric-models
include("../../test/simple_model.jl") #hide
````

![simple_model](./assets/simple_model.svg)

````@example 1-parametric-models
model
````

Build a basic ConstraintTree model without parameters

````@example 1-parametric-models
m = COBREXA.flux_balance_constraints(model)
````

Solve normally

````@example 1-parametric-models
base_model =
    COBREXA.optimized_values(m; optimizer = Tulip.Optimizer, objective = m.objective.value)
base_model.fluxes
````

## Add parameters to the model

Make bound of r2 and mass balance of m3 parameters

````@example 1-parametric-models
Symbolics.@variables r2bound m3bound

m.fluxes.r2 =
    ConstraintTrees.Constraint(m.fluxes.r2.value, -2 * ParameterBetween(r2bound, 0))

m.flux_stoichiometry.m3 =
    ConstraintTrees.Constraint(m.flux_stoichiometry.m3.value, ParameterEqualTo(m3bound) / 2)
````

# add parametric constraints

````@example 1-parametric-models
Symbolics.@variables p[1:4]

m *=
    :linparam^ConstraintTrees.Constraint(
        value = p[1] * m.fluxes.r1.value + p[2] * m.fluxes.r2.value,
        bound = -ParameterBetween(p[3], 0),
    )
````

substitute params into model

````@example 1-parametric-models
parameter_substitutions = Dict(
    r2bound => 4.0,
    m3bound => 0.1, # lose some mass here
    p[1] => 1.0,
    p[2] => 1.0,
    p[3] => 4.0,
)

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
)
m_noparams.fluxes
````

## Change the parameters and re-solve

substitute params into model

````@example 1-parametric-models
parameter_substitutions[m3bound] = 0.0

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Tulip.Optimizer,
)
m_noparams.fluxes
````

## Quadratic parameters also work

````@example 1-parametric-models
Symbolics.@variables q[1:6]

m.objective = ConstraintTrees.Constraint(
    value = sum(
        rxn.value * rxn.value * qi for (qi, rxn) in zip(collect(q), values(m.fluxes))
    ),
    bound = nothing,
)

m *= :objective_bound^ConstraintTrees.Constraint(value = m.fluxes.r6.value, bound = 2.0)

parameter_substitutions = merge(parameter_substitutions, Dict(zip(q, fill(1.0, 6))))

m_noparams, _, _, _ = optimized_constraints_with_parameters(
    m,
    parameter_substitutions;
    objective = m.objective.value,
    optimizer = Clarabel.Optimizer,
    sense = Minimal,
)
m_noparams.fluxes




pqv17 = ParameterLinearValue([1, 2], [2, 1]) * ParameterLinearValue([2], [1]) #src
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

