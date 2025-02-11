
#=
Copyright (c) 2025, Heinrich-Heine University Duesseldorf
Copyright (c) 2025, University of Luxembourg

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
    module DifferentiableMetabolism

Differentiate parameterized genome-scale metabolic models. Currently,
`FastDifferentiation` parameters are used in conjunction with `COBREXA` models.
To ensure your derivatives are well defined, you must ensure that the models are
*pruned*, i.e. all active reactions, genes, etc. have nonzero values at the
optimum. See the documentation for model details and examples.
"""
module DifferentiableMetabolism

using DocStringExtensions

import AbstractFBCModels as A
import COBREXA as X
import JuMP as J
import ConstraintTrees as C
import LinearAlgebra as LA
import SparseArrays as SA
import FastDifferentiation as F

const Ex = F.Node
const LinearValueP = C.LinearValueT{Ex}
const QuadraticValueP = C.QuadraticValueT{Ex}
const BetweenP = C.BetweenT{Ex}
const EqualToP = C.EqualToT{Ex}

include("parameter_promotion.jl")
include("substitute.jl")
include("misc.jl")
include("get_constraints.jl")
include("solver.jl")
include("differentiate.jl")
include("prune.jl")
include("public_api.jl")

end
