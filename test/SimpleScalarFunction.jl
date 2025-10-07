using Test
@testset "Simple Scalar Function" begin
    using FixedPointAcceleration

    func(x) = cos.(x)
    Inputs = 1.1
    fp_simple = fixed_point(func, Inputs, Simple())
    @test fp_simple.Convergence_ < 1e-10
    fp_anderson = fixed_point(func, Inputs, Anderson())
    @test fp_anderson.Convergence_ < 1e-10
    fp_aitken = fixed_point(func, Inputs, Aitken())
    @test fp_aitken.Convergence_ < 1e-10
    fp_newton = fixed_point(func, Inputs, Newton())
    @test fp_newton.Convergence_ < 1e-10
    fp_VEA = fixed_point(func, Inputs, VEA())
    @test fp_VEA.Convergence_ < 1e-10
    fp_SEA = fixed_point(func, Inputs, SEA())
    @test fp_SEA.Convergence_ < 1e-10
    fp_MPE = fixed_point(func, Inputs, MPE())
    @test fp_MPE.Convergence_ < 1e-10
    fp_RRE = fixed_point(func, Inputs, RRE())
    @test fp_RRE.Convergence_ < 1e-10

    # Testing input/output matrix functionality
    a = fixed_point(func, fp_RRE.Inputs_, fp_RRE.Outputs_, RRE(); MaxIter=20)
    @test a.Convergence_ < 1e-8 # If the above line throws then that is a test failure.
    # Testing different input/output shapes
    b = fixed_point(
        func, fp_RRE.Inputs_[:, 1:1], fp_RRE.Outputs_[:, 1:1], RRE(); MaxIter=20
    )
    @test b.Convergence_ < 1e-8 # If the above line throws then that is a test failure.
end
