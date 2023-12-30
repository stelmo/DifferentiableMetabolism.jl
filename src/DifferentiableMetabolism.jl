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
include("parameter_linearvalue.jl")
include("parameter_isozyme.jl")
include("parameter_quadraticvalue.jl")

# ConstraintTrees and Symbolics
include("constraint_trees.jl")
include("symbolics.jl")

# the juice
include("constraints.jl")
include("solver.jl")
include("kkt.jl")

end
