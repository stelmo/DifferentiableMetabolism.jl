module DifferentiableMetabolism

using DocStringExtensions

import COBREXA
import JuMP
import ConstraintTrees
import Symbolics
import LinearAlgebra
import SparseArrays: sparse, sparsevec

# Define new parameter-based types
include("parameter_bound.jl")
include("parameter_isozyme.jl")
include("parameter_linearvalue.jl")
include("parameter_quadraticvalue.jl")
include("parameter_promotion.jl")

# ConstraintTrees and Symbolics
include("constraint_trees.jl")
include("symbolics.jl")

# the juice
include("get_constraints.jl")
include("solver.jl")
include("kkt.jl")

end
