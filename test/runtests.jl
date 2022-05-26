using DifferentiableMetabolism
using Test
using Tulip, Ipopt
using COBREXA
using SparseArrays, LinearAlgebra

# Testing infrastructure taken from COBREXA
TEST_TOLERANCE = 1e-6
TEST_TOLERANCE_RELAXED = 1e-3

print_timing(fn, t) = @info "$(fn) done in $(round(t; digits = 2))s"

# Helper functions for running tests en masse
function run_test_file(path...)
    fn = joinpath(path...)
    t = @elapsed include(fn)
    print_timing(fn, t)
end

run_test_file("static_data.jl")

@testset "DifferentiableMetabolism.jl" begin
    run_test_file("enzyme.jl")
    run_test_file("prune.jl")
    run_test_file("update.jl")
    run_test_file("differentiable_models/basic_gecko.jl")
    run_test_file("differentiable_models/basic_smoment.jl")
    run_test_file("differentiable_models/thermodynamic_smoment.jl")
    run_test_file("differentiable_models/thermodynamic_gecko.jl")
    run_test_file("differentiable_models/michaelis_menten_gecko.jl")
    run_test_file("analytic_derivatives_with_scaling.jl")
end
