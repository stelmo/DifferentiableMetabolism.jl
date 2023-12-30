
import DifferentiableMetabolism
using Test

const TEST_TOLERANCE = 1e-3

@testset "DifferentiableMetabolism tests" begin
    @testset "Parametric models" begin
        include("../docs/src/1-parametric-models.jl")
    end

    @testset "Differentiating enzyme constrained models" begin
        include("../docs/src/2-differentiate-enzyme-model.jl")
    end

    @testset "Differentiating nonlinear models" begin
        include("../docs/src/3-differentiate-nonlinear-model.jl")
    end
end
