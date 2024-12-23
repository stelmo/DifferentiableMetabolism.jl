#=
Copyright (c) 2024, Heinrich-Heine University Duesseldorf
Copyright (c) 2024, University of Luxembourg

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

import AbstractFBCModels
import COBREXA
import JuMP
import ConstraintTrees
import ConstraintTrees: LinearValue, LinearValueT, QuadraticValue, QuadraticValueT, Between, BetweenT, EqualTo, EqualToT
import LinearAlgebra
import SparseArrays
import FastDifferentiation

const Expression = FastDifferentiation.Node

# Define new parameter-based types
include("parameter_promotion.jl")

# ConstraintTrees and Expression handling utils
include("substitute.jl")

# the juice
include("get_constraints.jl")
include("solver.jl")
# include("differentiate.jl")
# include("differentiate_efm.jl")
# include("prune.jl")

end
