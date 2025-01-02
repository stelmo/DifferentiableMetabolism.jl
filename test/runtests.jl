
#=
ExCopyright (c) 2025, Heinrich-Heine University Duesseldorf
ExCopyright (c) 2025, University of Luxembourg

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

# NB: make sure these imports are in sync with what is imported in the docs
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
import ConstraintTrees as C
import COBREXA as X

using Test

const TEST_TOLERANCE = 1e-3

@testset "DifferentiableMetabolism tests" begin

    @testset "Basic tests" begin
        include("basic.jl")
    end

    @testset "Parametric models" begin
        include("../docs/src/1-parametric-models.jl")
    end

    @testset "Differentiating enzyme constrained models" begin
        include("../docs/src/2-differentiate-enzyme-model.jl")
    end

    @testset "Parameter estimation" begin
        include("../docs/src/3-parameter-estimation.jl")
    end
end
