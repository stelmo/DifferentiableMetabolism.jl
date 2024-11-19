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
# cost matrix of the EFMs
parameters = 1
w_jk = molar_mass(gene_j)/pool_k
V_ij = flux_j_in_efm_i
kcat_j = kcat(reaction_j)
d_ik = sum()


f = sum(lambda)
g = D*lambda - ones(size(D,1))


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
    EFMs,
    θ,
    optimizer
)
    n_vars = size(calc_D(θ,EFMs),2)
    efm_opt = JuMP.Model(optimizer)
    @variable(efm_opt, z[1:n_vars])
    @constraint(efm_opt, eq, float.(calc_D(θ,EFMs))*z ==[1 ; 1])
    @objective(efm_opt, Max, sum(z))
    optimize!(efm_opt)

    x = value.(efm_opt[:z]) 
    ν = dual.(efm_opt[:eq])

    # define L, the gradient of the Lagrangian 
    L(x,ν,θ) = [
        ones(n_vars) + calc_D(θ,EFMs)'*ν
        calc_D(θ,EFMs)*x - ones(n_vars)
    ]
    # differentiate L wrt x,ν, the variables
    dL_vars(x,ν,θ) = [
        spzeros(n_vars,n_vars) calc_D(θ,EFMs)'
        calc_D(θ,EFMs)   spzeros(n_vars,n_vars)
    ] 
    # differentiate L wrt θ
    dL_params(x,ν,θ) = ForwardDiff.jacobian(θ -> L(x, ν, θ), θ)
    # solve for d_vars/d_params 
    dx = -Array(dL_vars(x,ν,θ))\dL_params(x,ν,θ)
    
    # note: dx[[3,4],:] gives the derivatives of the dual variables ν 
    return dx[[1,2],:]
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
    θ,
    gene_product_mass_group,
    gene_product_mass_group_bound,
    gene_product_molar_mass,
    rid_enzyme,
)
    D = Matrix{Real}(undef,length(gene_product_mass_group_bound),length(EFMs))
    for (i,(enzyme_group,enzyme_bound)) in enumerate(gene_product_mass_group_bound)
        for (j,efm) in enumerate(EFMs)
            D[i,j] = 0
            for (k,rid) in enumerate(collect(keys(rid_enzyme)))
                !any(((g,c),) -> gene_product_mass_group[g] == enzyme_group, rid_enzyme[rid].gene_product_count) && continue 

                D[i,j] += sum(
                    [
                        efm[rid] * gene_product_molar_mass[g] * c / (enzyme_bound * θ[k])
                        for (g,c) in rid_enzyme[rid].gene_product_count if gene_product_mass_group[g] == enzyme_group
                    ]
                )
            end
        end
    end
    return D
end

export cost_matrix