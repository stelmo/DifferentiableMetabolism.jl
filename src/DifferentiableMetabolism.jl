module DifferentiableMetabolism

using COBREXA, JuMP
using ForwardDiff
using LinearAlgebra, SparseArrays, RowEchelon

include("utils.jl")
include("gecko.jl")
include("smoment.jl")
include("thermodynamic_smoment.jl")
include("diff.jl")

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end


end
