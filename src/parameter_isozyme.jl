
#=
Copyright (c) 2023, Heinrich-Heine University Duesseldorf
Copyright (c) 2023, University of Luxembourg

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

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
