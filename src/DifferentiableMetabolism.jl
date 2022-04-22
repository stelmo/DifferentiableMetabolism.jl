module DifferentiableMetabolism

using COBREXA, JuMP
using ForwardDiff
using LinearAlgebra, SparseArrays, RowEchelon
using DocStringExtensions

include("DifferentiableModel.jl")
include("Enzyme.jl")
include("differentiate.jl")
include("scale.jl")
include("utils.jl")
include("gecko.jl")
include("gecko_analytic_derivatives.jl")
include("smoment.jl")
include("thermodynamic_smoment.jl")
include("thermodynamic_gecko.jl")
include("michaelis_menten_gecko.jl")

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end
