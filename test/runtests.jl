using FixedPointAccelerationNext
using Test

using TestItems

@testitem "quality assurance" begin
    include("quality_assurance.jl")
end

@testitem "Simple scalar function" begin
    include("simple_scalar_function.jl")
end
@testitem "simple vector function" begin
    include("simple_vector_function.jl")
end
@testitem "autodiff" begin
    include("autodiff.jl")
end
@testitem "acceleration_tests" begin
    include("acceleration_tests.jl")
end
@testitem "complex" begin
    include("complex.jl")
end
@testitem "allocations" begin
    include("test_allocations.jl")
end

# @testitem "consumption smoothing" begin
#     include("ConsumptionSmoothingProblem.jl")
# end

# @testitem "consumption smoothing old" begin
#     include("ConsumptionSmoothingProblem_old.jl")
# end

using TestItemRunner

@run_package_tests verbose=true

# @testset "test_side_products" include("test_side_products.jl")
