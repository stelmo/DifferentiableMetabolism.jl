# # Differentiating the EFMs/OFMs of an optimal solution

# The optimal flux distribution of any metabolic model can be written as a 
# weighted sum of the EFMs of that model. We are interested in calculating 
# the sensitivity of these weightings to the model parameters. 

using DifferentiableMetabolism

using ElementaryFluxModes
using FastDifferentiation
using COBREXA
using Tulip 
using ConstraintTrees
using LinearAlgebra

# ## Build a simple enzyme constrained model

# The code used to construct the model is located in `test/simple_model.jl`, but
# it is not shown here for brevity.

include("../../test/efm_simple_model.jl"); #hide

model

parameter_values = Dict{Symbol, Float64}()
reaction_isozymes = Dict{String,Dict{String,ParameterIsozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
float_reaction_isozymes = Dict{String,Dict{String,COBREXA.Isozyme}}() #src
for rid in AbstractFBCModels.reactions(model)
    grrs = AbstractFBCModels.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(rid_kcat, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)

        kcat = rid_kcat[rid]
        parameter_values[Symbol(rid)] = kcat

        d = get!(reaction_isozymes, rid, Dict{String,ParameterIsozyme}())
        d["isozyme_$i"] = ParameterIsozyme(
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = FastDifferentiation.Node(Symbol(rid)),
            kcat_reverse = 0.0,
        )
        d2 = get!(float_reaction_isozymes, rid, Dict{String,COBREXA.Isozyme}()) #src
        d2["isozyme_$i"] = COBREXA.Isozyme( #src
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), #src
            kcat_forward = kcat, #src
            kcat_reverse = 0.0, #src
        ) #src
    end
end


km = build_kinetic_model(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
)

ec_solution, _, _, _ = optimized_constraints_with_parameters(
    km,
    parameter_values;
    objective = km.objective.value,
    optimizer = Tulip.Optimizer,
    modifications = [COBREXA.set_optimizer_attribute("IPM_IterationsLimit", 10_000)],
)

ec_solution.fluxes

ec_solution_fba = enzyme_constrained_flux_balance_analysis( #src
    model; #src
    reaction_isozymes = float_reaction_isozymes, #src
    gene_product_molar_masses, #src
    capacity, #src
    optimizer = Tulip.Optimizer, #src
) #src

@test isapprox(ec_solution.objective, ec_solution_fba.objective; atol = TEST_TOLERANCE) #src

@test any(isapprox.(values(ec_solution.gene_product_amounts), 0, atol=1e-8)) #src

# We have a solution that uses every reaction, and the enzyme capacities are both full.
# Therefore, we may calculate the EFMs of this solution and direactly differentiate
# them, with no pruning required.

# ## Calculate EFMs of the optimal solution 

# We need to input the stoichiometric matrix `N`, and a nullspace `K`, into 
# ElementaryFluxModes.jl

N = AbstractFBCModels.stoichiometry(model)

K = rational_nullspace(Matrix(N))[1]

# Permute the rows of `K` to be in the form `[I;K*]`

order = Int64[]
rows_done = 1
for (i, row) in enumerate(eachrow(K))
    global rows_done
    rows_done > size(K, 2) && break
    if row == Matrix(I(size(K, 2)))[rows_done, :]
        push!(order, i)
        rows_done += 1
    end
end
append!(order, [i for i = 1:size(K, 1) if i ∉ order])
K = K[order, :]

# The reaction order of `N` must match that of `K`
N = N[:, order]

# Run the double description algorithm 

R = DDBinary(N, K)
E = Matrix(undef, size(R, 1), size(R, 2))
for (i, r) in enumerate(eachcol(R))
    non_zero = findall(x -> x != 0, r)
    flux_ns = rational_nullspace(Matrix(N[:, non_zero]); tol = 1e-14)[1]
    mode = zeros(size(R, 1))
    for (j, x) in zip(non_zero, flux_ns)
        mode[j] = abs(x) < 1e-14 ? 0 : x
    end
    E[:, i] = mode
end
E = E[invperm(order), :] # put the entries into original reaction order

# Make a dictionary out of the EFM result 

EFM_dict = Dict(AbstractFBCModels.reactions(model) .=> eachrow(E))
EFMs = [
    Dict(k => v[1] / EFM_dict["r6"][1] for (k, v) in EFM_dict),
    Dict(k => v[2] / EFM_dict["r6"][2] for (k, v) in EFM_dict)
]

# ## Differentiate the EFMs 

# We have calculated the EFMs, and now wish to differentiate their weightings 
# with respect to the model parameters.

parameters = FastDifferentiation.Node.(collect(keys(parameter_values)))
p_vals = collect(values(parameter_values))
rid_pid = Dict(rid => [iso.kcat_forward for (k, iso) in v][1] for (rid, v) in reaction_isozymes)
rid_gcounts = Dict(rid => [v.gene_product_stoichiometry for (k, v) in d][1] for (rid, d) in reaction_isozymes)
sens_efm = differentiate_efm(EFMs, parameters, rid_pid, parameter_values, rid_gcounts, capacity, gene_product_molar_masses, Tulip.Optimizer)
