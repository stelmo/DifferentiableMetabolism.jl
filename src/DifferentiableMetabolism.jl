module DifferentiableMetabolism

using COBREXA, JuMP
using ForwardDiff, Symbolics
using LinearAlgebra, SparseArrays, RowEchelon
using DocStringExtensions

include("DifferentiableModel.jl")
include("Enzyme.jl")
include("differentiate.jl")
include("scale.jl")
include("utils.jl")
include("update.jl")
include("gecko.jl")
include("smoment.jl")
include("thermodynamic_smoment.jl")
include("thermodynamic_gecko.jl")
include("michaelis_menten_gecko.jl")
include("analytic_derivatives.jl")
include("inplace_differentiate.jl")
include("inplace_derivatives.jl")

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end
