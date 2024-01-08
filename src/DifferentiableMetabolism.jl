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

module DifferentiableMetabolism

using DocStringExtensions

import COBREXA
import JuMP
import ConstraintTrees
import Symbolics
import LinearAlgebra: :(\), rank, qr
import SparseArrays: sparse, sparsevec, findnz, dropzeros, dropzeros!

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
include("differentiate.jl")

end
