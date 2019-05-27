using FixedPointAcceleration
using Test

# Run tests
println("Test putting together iterates with jumps")
@time @test include("TestPuttingInputsAndOutputsTogether.jl")
println("Test simple Scalar function")
@time @test include("SimpleScalarFunction.jl")
println("Test simple vector function")
@time @test include("SimpleVectorFunction.jl")
println("Test Chaining together calls with cos functions.")
@time @test include("Chaining.jl")
println("Testing Bounds Error")
@time @test include("BoundsError.jl")
println("Testing ReplaceInvalidsTest")
@time @test include("ReplaceInvalidsTesting.jl")
println("CommonIncrementsTest")
@time @test include("CommonIncrementsTest.jl")
# Note consumption smoothing problem is not tested for time and dependency reasons.
# Note that the Autodiff extention is not included for the same reason.
