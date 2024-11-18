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
parameters = 
w_jk = molar_mass(gene_j)/pool_k
V_ij = flux_j_in_efm_i
kcat_j = kcat(reaction_j)
d_ik = sum()


f = sum(lambda)
g = D*lambda - ones(size(D,1))