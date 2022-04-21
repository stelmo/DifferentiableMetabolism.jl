using DifferentiableMetabolism
using Test
using Tulip, Ipopt
using COBREXA
using JSON

# Testing infrastructure taken from COBREXA
TEST_TOLERANCE = 1e-6

print_timing(fn, t) = @info "$(fn) done in $(round(t; digits = 2))s"

# Helper functions for running tests en masse
function run_test_file(path...)
    fn = joinpath(path...)
    t = @elapsed include(fn)
    print_timing(fn, t)
end

run_test_file("static_data.jl")

@testset "DifferentiableMetabolism.jl" begin
    run_test_file("gecko.jl")
    # run_test_file("smoment.jl")
    # run_test_file("thermodynamic_smoment.jl")
end
