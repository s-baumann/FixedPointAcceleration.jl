using FixedPointAcceleration
using Test

# Run tests
println("Test putting together iterates with jumps")
@time @test include("TestPuttingInputsAndOutputsTogether.jl")
println("Test simple Scalar function")
@time @test include("SimpleScalarFunction.jl")
println("Test simple vector function")
@time @test include("SimpleVectorFunction.jl")
