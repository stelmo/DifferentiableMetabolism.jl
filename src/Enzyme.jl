"""
    $(TYPEDEF)

A struct used to store information about an enzyme. A simplified version of
[`COBREXA.Isozyme`](@ref), which stores only the turnover number (kcat) that
should be used by the differentiable model, as well as information about the
protein complex (stoichiometry and molar mass).

# Fields
$(TYPEDFIELDS)
"""
mutable struct Enzyme
    kcat::Float64
    gene_product_count::Dict{String,Int}
    gene_product_mass::Dict{String,Float64}
end

"""
    $(TYPEDSIGNATURES)

Convert a [`COBREXA.Isozyme`](@ref) to an [`Enzyme`](@ref) using the turnover
number in `direction` (either `:forward` or `:reverse`), as well as
`gene_product_mass_lookup`, which is either a function or a dictionary mapping a
gene product id to its molar mass.
"""
function isozyme_to_enzyme(
    isozyme::Isozyme,
    gene_product_mass_lookup::Union{Function,Dict{String,Float64}};
    direction = :forward,
)
    _gpmlu =
        gene_product_mass_lookup isa Function ? gene_product_mass_lookup :
        (gid -> gene_product_mass_lookup[gid])
    
    Enzyme(
        direction == :forward ? isozyme.kcat_forward : isozyme.kcat_reverse,
        isozyme.gene_product_count,
        Dict(gid => _gpmlu(gid) for gid in keys(isozyme.gene_product_count)),
    )
end
