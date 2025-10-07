module test_simple_scalar

using Test
@testset "Simple Scalar Function" begin
    using FixedPointAcceleration

    func(x) = cos.(x)
    Inputs = 1.1
    fp_simple = fixed_point(func, Inputs, Simple())
    @test fp_simple.convergence < 1e-10
    fp_anderson = fixed_point(func, Inputs, Anderson())
    @test fp_anderson.convergence < 1e-10
    fp_aitken = fixed_point(func, Inputs, Aitken())
    @test fp_aitken.convergence < 1e-10
    fp_newton = fixed_point(func, Inputs, Newton())
    @test fp_newton.convergence < 1e-10
    fp_VEA = fixed_point(func, Inputs, VEA())
    @test fp_VEA.convergence < 1e-10
    fp_SEA = fixed_point(func, Inputs, SEA())
    @test fp_SEA.convergence < 1e-10
    fp_MPE = fixed_point(func, Inputs, MPE())
    @test fp_MPE.convergence < 1e-10
    fp_RRE = fixed_point(func, Inputs, RRE())
    @test fp_RRE.convergence < 1e-10

    # Testing input/output matrix functionality
    opts_20 = FixedPointOptions(max_iterations=20)
    a = fixed_point(func, fp_RRE.inputs, fp_RRE.outputs, RRE(), opts_20)
    @test a.convergence < 1e-8 # If the above line throws then that is a test failure.
    # Testing different input/output shapes
    b = fixed_point(func, fp_RRE.inputs[:, 1:1], fp_RRE.outputs[:, 1:1], RRE(), opts_20)
    @test b.convergence < 1e-8 # If the above line throws then that is a test failure.
end
end
