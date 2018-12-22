using FixedPoint
using Base.Test

# Run tests

tic()
println("Test simple vector function")
@time @test include("SecondDerivativeTest.jl")
toc()
