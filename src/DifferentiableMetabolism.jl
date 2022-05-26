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
include("differentiable_models/gecko.jl")
include("differentiable_models/smoment.jl")
include("analytic_derivatives/analytic_derivatives.jl")
include("analytic_derivatives/analytic_derivatives_with_scaling.jl")

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end
