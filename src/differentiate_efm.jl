"""
$(TYPEDSIGNATURES)

Differentiate a set of elementary flux modes `EFMs` with respect to `parameters` which take on 
values `parameter_values` in the optimal solution with respect to the optimal usage of these EFMs 
in a solution. 

Calculating the proportion of each EFM used in the optimal solution can be turned into the 
following LP:
    maximise     ∑xᵢ 
    subject to   Dx = 1
where D is the cost matrix associated to each EFM in each enzyme pool.   
We then differentiate Lagrangian of this LP to calculate the differential of x by parameters.

"""
# # cost matrix of the EFMs
# parameters = 1
# w_jk = molar_mass(gene_j)/pool_k
# V_ij = flux_j_in_efm_i
# kcat_j = kcat(reaction_j)
# d_ik = sum()


# f = sum(lambda)
# g = D*lambda - ones(size(D,1))


"""
$(TYPEDSIGNATURES)

Differentiate the usage of EFMs with respect to changes in the model parameters using implicit 
differentiation of the Lagrangian.
Calculating the proportion of each EFM used in the optimal solution can be turned into the 
following LP:
    maximise     ∑xᵢ 
    subject to   Dx = 1
where D is the cost matrix associated to each EFM in each enzyme pool.   
We then differentiate Lagrangian of this LP to calculate the differential of x by parameters.

Variables:
- 'D_matrix': the matrix of cost vectors, and must be inputted as a function of the parameters 
- 'parameters': vector of the model parameters
"""
function differentiate_efm(
    EFMs::Vector{Dict{String,Float64}},
    parameters,
    rid_pid,
    parameter_values,
    rid_gcounts,
    capacity::Vector{Tuple{String,Vector{String},Float64}},
    gene_product_molar_masses::Dict{String,Float64},
    optimizer
)
    D_eval = float.(cost_matrix(
        EFMs,
        rid_pid,
        rid_gcounts,
        capacity,
        gene_product_molar_masses;
        evaluate=true,
        parameter_values
    ))
    n_vars = size(D_eval, 2)

    efm_opt = JuMP.Model(optimizer)
    JuMP.@variable(efm_opt, z[1:n_vars])
    JuMP.@constraint(efm_opt, eq, D_eval * z == [1; 1])


    JuMP.@objective(efm_opt, Max, sum(z))
    JuMP.optimize!(efm_opt)

    x = JuMP.value.(efm_opt[:z])
    ν = JuMP.dual.(efm_opt[:eq])
    D = FastDifferentiation.Node.(
        cost_matrix(
            EFMs,
            rid_pid,
            rid_gcounts,
            capacity,
            gene_product_molar_masses,
        )
    )

    # define L, the gradient of the Lagrangian 
    L(x, ν, parameters) = [
        ones(n_vars) + D' * ν
        D * x - ones(n_vars)
    ]
    # differentiate L wrt x,ν, the variables
    dl_vars = [
        SparseArrays.spzeros(n_vars, n_vars) D'
        D SparseArrays.spzeros(n_vars, n_vars)
    ]

    # differentiate L wrt parameters
    dL_params(x, ν, parameters) = FastDifferentiation.jacobian(L(x, ν, parameters), parameters)
    # substitute parameter values: 
    dL_params_eval = make_function(dL_params(x, ν, parameters), parameters)

    dx = -Array(dl_vars) \ dL_params_eval(param_vals)

    # note: dx[[3,4],:] gives the derivatives of the dual variables ν 
    return dx[[1, 2], :]
end

export differentiate_efm

"""
$(TYPEDSIGNATURES)

Calculate a matrix of the cost vectors of each EFM to each constraint.
Entry (i,j) gives the total cost in constraint i to produce one unit objective flux through EFM j.
Cost is calculated as ∑w(i)V(j)/kcat, where the variables are:
- 'w(i)': fraction of the ith enzyme pool that one mole of the enzyme uses up
- 'V(j)': flux through the reaction in EFM j 
- 'kcat': turnover number of the enzyme.

Inputted function variables are:
- 'efms': list of the fluxes through the EFMs, each given as a dictionary of reaction_id => [flux efm1, flux efm2, ...]
- 'parameters': parameters of the LP, the turnover numbers 
"""

function cost_matrix(
    EFMs::Vector{Dict{String,Float64}},
    rid_pid,
    rid_gcounts,
    capacity::Vector{Tuple{String,Vector{String},Float64}},
    gene_product_molar_masses::Dict{String,Float64};
    evaluate=false,
    parameter_values=nothing
)
    D = Matrix(undef, length(capacity), length(EFMs))
    for (i, (enzyme_group, enzymes, enzyme_bound)) in enumerate(capacity)
        for (j, efm) in enumerate(EFMs)
            D[i, j] = 0
            for (rid, gcount) in rid_gcounts
                pid = rid_pid[rid]

                # if genes not in pool i, skip 
                all(((g, c),) -> g ∉ enzymes, gcount) && continue

                # otherwise, add the cost of this enzyme to the ith pool from the jth efm
                if !evaluate
                    D[i, j] += sum(
                        [
                        efm[rid] * gene_product_molar_masses[g] * c / (enzyme_bound * pid)
                        for (g, c) in gcount if g ∈ enzymes
                    ]
                    )
                else
                    D[i, j] += Float64(sum(
                        [
                        efm[rid] * gene_product_molar_masses[g] * c / (enzyme_bound * parameter_values[Symbol(pid)])
                        for (g, c) in gcount if g ∈ enzymes
                    ]
                    ))
                end
            end
        end
    end
    return D
end

export cost_matrix