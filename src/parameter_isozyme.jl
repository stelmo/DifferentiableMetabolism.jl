
@kwdef mutable struct ParameterIsozyme{T} <: X.Isozyme
    gene_product_stoichiometry::Dict{String,T} # stoichiometry could also be parameters, but not required
    kcat_forward::X.Maybe{S.Num} = nothing
    kcat_backward::X.Maybe{S.Num} = nothing
end

export ParameterIsozyme
