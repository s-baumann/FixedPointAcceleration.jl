using FixedPointAccelerationNext
using Test

using TestItems

@testitem "quality assurance" begin
    include("Aqua.jl")
end

@testitem "SimpleScalarFunction" begin
    include("SimpleScalarFunction.jl")
end
@testitem "SimpleVectorFunction" begin
    include("SimpleVectorFunction.jl")
end
@testitem "autodiff" begin
    include("Autodiff.jl")
end
@testitem "new_core" begin
    include("test_new_core.jl")
end
@testitem "new_core_complex" begin
    include("test_new_core_complex.jl")
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

# @testset "ComplexNumberTests" include("ComplexNumberTests.jl")
# @testset "parity_old_new" include("test_parity_old_new.jl")
# @testset "test_side_products" include("test_side_products.jl")
# @testset "input_output" include("TestPuttingInputsAndOutputsTogether.jl")
# @testset "CommonIncrements" include("CommonIncrementsTest.jl")
# @testset "Chaining" include("Chaining.jl")
# @testset "BoundsError" include("BoundsError.jl")
# @testset "ReplaceInvalidsTesting" include("ReplaceInvalidsTesting.jl")
