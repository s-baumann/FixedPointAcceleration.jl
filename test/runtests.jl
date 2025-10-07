using FixedPointAcceleration
using Test

# Run tests
@testset "quality assurance" include("Aqua.jl")
@testset "input_output" include("TestPuttingInputsAndOutputsTogether.jl")
@testset "CommonIncrements" include("CommonIncrementsTest.jl")
@testset "SimpleScalarFunction" include("SimpleScalarFunction.jl")
@testset "SimpleVectorFunction" include("SimpleVectorFunction.jl")
@testset "Chaining" include("Chaining.jl")
@testset "BoundsError" include("BoundsError.jl")
@testset "ReplaceInvalidsTesting" include("ReplaceInvalidsTesting.jl")
@testset "ComplexNumberTests" include("ComplexNumberTests.jl")
# @testset "test_side_products" include("test_side_products.jl")
@testset "consumption smoothing" include("ConsumptionSmoothingProblem.jl")
