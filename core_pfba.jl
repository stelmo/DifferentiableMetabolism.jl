

# add small quadratic weight
m.objective = ConstraintTrees.Constraint(
    value = 1e-5 * (
        sum(x.value * x.value for x in values(m.fluxes)) +
        sum(x.value * x.value for x in values(m.:gene_product_amounts))
    ) - m.fluxes.BIOMASS_Ecoli_core_w_GAM.value,
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

ec_solution.fluxes.BIOMASS_Ecoli_core_w_GAM
ec_solution.fluxes
ec_solution.enzymes

sens = differentiate(
    m,
    m.objective.value,
    x_vals,
    eq_dual_vals,
    ineq_dual_vals,
    parameter_values,
    [capacitylimitation; kcats],
)

objective = m.objective.value
parameters = [capacitylimitation; kcats]
