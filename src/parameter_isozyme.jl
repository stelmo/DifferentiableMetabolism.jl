
"""
$(TYPEDEF)

Subtype of [`COBREXA.Isozyme`](@ref) which includes parameters in
`gene_product_stoichiometry` (optionally), the `kcat_forward`, and the
`kcat_backward`. If the reaction does not have a turnover number, `nothing` can
be used. 

# Fields
$(TYPEDFIELDS)
"""
@kwdef mutable struct ParameterIsozyme{T} <: COBREXA.Isozyme
    gene_product_stoichiometry::Dict{String,T}
    kcat_forward::COBREXA.Maybe{Symbolics.Num} = nothing
    kcat_backward::COBREXA.Maybe{Symbolics.Num} = nothing
end

export ParameterIsozyme
