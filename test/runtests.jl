using FixedPointAcceleration
using Test

@testset "FixedPointAcceleration.jl" begin
    @testset "quality assurance" include("Aqua.jl")

    @testset "SimpleScalarFunction" include("SimpleScalarFunction.jl")
    @testset "SimpleVectorFunction" include("SimpleVectorFunction.jl")

    @testset "new_core" include("test_new_core.jl")
    @testset "new_core_complex" include("test_new_core_complex.jl")
    @testset "allocations" include("test_allocations.jl")
end

# @testset "ComplexNumberTests" include("ComplexNumberTests.jl")
# @testset "parity_old_new" include("test_parity_old_new.jl")
# @testset "test_side_products" include("test_side_products.jl")
# @testset "consumption smoothing" include("ConsumptionSmoothingProblem.jl")
# @testset "input_output" include("TestPuttingInputsAndOutputsTogether.jl")
# @testset "CommonIncrements" include("CommonIncrementsTest.jl")
# @testset "Chaining" include("Chaining.jl")
# @testset "BoundsError" include("BoundsError.jl")
# @testset "ReplaceInvalidsTesting" include("ReplaceInvalidsTesting.jl")
