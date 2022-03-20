using FixedPointAcceleration
using Test

# Run tests
println("Test putting together iterates with jumps")
include("TestPuttingInputsAndOutputsTogether.jl")
println("Test simple Scalar function")
include("SimpleScalarFunction.jl")
println("Test simple vector function")
include("SimpleVectorFunction.jl")
println("Test Chaining together calls with cos functions.")
include("Chaining.jl")
println("Testing Bounds Error")
include("BoundsError.jl")
println("Testing ReplaceInvalidsTest")
include("ReplaceInvalidsTesting.jl")
#println("CommonIncrementsTest")
#@time @test include("CommonIncrementsTest.jl")
#println("test_side_products")
#@time @test include("test_side_products.jl")

# Note consumption smoothing problem is not tested for time and dependency reasons.
# Note that the Autodiff extention is not included for the same reason.
# And test_side_products for the same reason.
