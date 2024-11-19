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
We then differentiate Lagrangian of this LP to calculate the differential of x by θ.

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
We then differentiate Lagrangian of this LP to calculate the differential of x by θ.

Variables:
- 'D_matrix': the matrix of cost vectors, and must be inputted as a function of the parameters 
- 'θ': vector of the model parameters
"""
function differentiate_efm(
    EFMs::Vector{Dict{String, Float64}},
    θ::Vector{Symbol},
    reaction_parameter_isozymes::Dict{String,Dict{String,ParameterIsozyme}},
    capacity::Vector{Tuple{String,Vector{String},Float64}},
    gene_product_molar_masses::Dict{String,Float64},
    optimizer
)

    D(θ) = cost_matrix(
        EFMs,
        reaction_parameter_isozymes,
        capacity,
        gene_product_molar_masses,
    )
    n_vars = size(D(θ), 2)
    efm_opt = JuMP.Model(optimizer)
    JuMP.@variable(efm_opt, z[1:n_vars])
    JuMP.@constraint(efm_opt, eq, float.(D(θ)) * z == [1; 1])
    JuMP.@objective(efm_opt, Max, sum(z))
    optimize!(efm_opt)

    x = value.(efm_opt[:z])
    ν = dual.(efm_opt[:eq])



    # define L, the gradient of the Lagrangian 
    L(x, ν, θ) = [
        ones(n_vars) + D(θ)' * ν
        D(θ) * x - ones(n_vars)
    ]
    # differentiate L wrt x,ν, the variables
    dL_vars(x, ν, θ) = [
        spzeros(n_vars, n_vars) D(θ)'
        D(θ) spzeros(n_vars, n_vars)
    ]
    # differentiate L wrt θ
    dL_params(x, ν, θ) = FastDifferentiation.jacobian(θ -> L(x, ν, θ), θ)
    # solve for d_vars/d_params 
    dx = -Array(dL_vars(x, ν, θ)) \ dL_params(x, ν, θ)

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
- 'θ': parameters of the LP, the turnover numbers 
"""
function cost_matrix(
    EFMs::Vector{Dict{String,Float64}},
    reaction_parameter_isozymes::Dict{String,Dict{String,ParameterIsozyme}},
    capacity::Vector{Tuple{String,Vector{String},Float64}},
    gene_product_molar_masses::Dict{String,Float64},
)
    parameter_values = Dict{Symbol,Float64}()
    for (r, iso) in reaction_parameter_isozymes
        for (k, v) in iso
            parameter_values[Symbol(v.kcat_forward)] = reaction_isozymes[r][k].kcat_forward
        end
    end
    rid_gcounts = Dict(rid => [v.gene_product_stoichiometry for (k, v) in d][1] for (rid, d) in reaction_parameter_isozymes)
    rid_pid = Dict(rid => [Symbol(iso.kcat_forward) for (k, iso) in v][1] for (rid, v) in reaction_parameter_isozymes)

    D = Matrix{Real}(undef, length(capacity), length(EFMs))
    for (i, (enzyme_group, enzymes, enzyme_bound)) in enumerate(capacity)
        for (j, efm) in enumerate(EFMs)
            D[i, j] = 0
            for (rid, gcount) in rid_gcounts
                pid = rid_pid[rid]

                # if genes not in pool i, skip 
                all(((g, c),) -> g ∉ enzymes, gcount) && continue

                # otherwise, add the cost of this enzyme to the ith pool from the jth efm
                D[i, j] += sum(
                    [
                    efm[rid] * gene_product_molar_masses[g] * c / (enzyme_bound * parameter_values[pid])
                    for (g, c) in gcount if g ∈ enzymes
                ]
                )
            end
        end
    end
    return D
end

export cost_matrix