module test_complex

using Test
@testset "Complex Number Functions" begin
    using FixedPointAccelerationNext.OldImplementation:
        fixed_point, Simple, Anderson, Aitken, Newton, MPE, RRE, VEA, SEA, FixedPointOptions
    # Core already loaded at top-level in runtests; avoid duplicate using

    # Define option configurations used throughout tests
    opts_100 = FixedPointOptions(max_iterations=100)
    opts_10 = FixedPointOptions(max_iterations=10)
    opts_silent = FixedPointOptions(quiet_errors=true)

    @testset "Simple Complex Scalar" begin
        # Test a simple complex function with a known fixed point
        func(x) = [0.5 * x[1] + 0.25im]
        Inputs = [1.0 + 0.5im]

        # Test with different algorithms
        fp_simple = fixed_point(func, Inputs, Simple())
        @test fp_simple.convergence < 1e-10
        @test abs(fp_simple.fixed_point[1] - (0.0 + 0.5im)) < 1e-10

        fp_anderson = fixed_point(func, Inputs, Anderson())
        @test fp_anderson.convergence < 1e-10
        @test abs(fp_anderson.fixed_point[1] - (0.0 + 0.5im)) < 1e-10

        fp_aitken = fixed_point(func, Inputs, Aitken())
        @test fp_aitken.convergence < 1e-10
        @test abs(fp_aitken.fixed_point[1] - (0.0 + 0.5im)) < 1e-10

        fp_newton = fixed_point(func, Inputs, Newton())
        @test fp_newton.convergence < 1e-10
        @test abs(fp_newton.fixed_point[1] - (0.0 + 0.5im)) < 1e-10

        fp_VEA = fixed_point(func, Inputs, VEA())
        @test fp_VEA.convergence < 1e-10
        @test abs(fp_VEA.fixed_point[1] - (0.0 + 0.5im)) < 1e-10

        fp_SEA = fixed_point(func, Inputs, SEA())
        @test fp_SEA.convergence < 1e-10
        @test abs(fp_SEA.fixed_point[1] - (0.0 + 0.5im)) < 1e-10

        fp_MPE = fixed_point(func, Inputs, MPE())
        @test fp_MPE.convergence < 1e-10
        @test abs(fp_MPE.fixed_point[1] - (0.0 + 0.5im)) < 1e-10

        fp_RRE = fixed_point(func, Inputs, RRE())
        @test fp_RRE.convergence < 1e-10
        @test abs(fp_RRE.fixed_point[1] - (0.0 + 0.5im)) < 1e-10
    end

    @testset "Complex Vector Functions" begin
        # Test a 2D complex function
        func(x) = [0.5 * x[1] + 0.1im, 0.3 * x[2] - 0.2im]
        Inputs = [1.0 + 0.2im, 0.5 - 0.3im]

        fp_simple = fixed_point(func, Inputs, Simple())
        @test fp_simple.convergence < 1e-10

        fp_anderson = fixed_point(func, Inputs, Anderson())
        @test fp_anderson.convergence < 1e-10

        # Test that the fixed point satisfies f(x) = x
        if !ismissing(fp_anderson.fixed_point)
            residual = func(fp_anderson.fixed_point) .- fp_anderson.fixed_point
            @test maximum(abs.(residual)) < 1e-10
        end
    end

    @testset "Complex Scalar Input" begin
        # Test with single complex number input (not vector)
        func(x) = [0.8 * x[1] + 0.1im]
        result = fixed_point(func, 1.0 + 0.5im, Anderson())
        @test result.convergence < 1e-10
        @test !ismissing(result.fixed_point)
    end

    @testset "Complex with Real Output" begin
        # Test a function that takes complex input but produces complex output (but real part dominates)
        func(x) = [0.8 * x[1] + 0.1]  # Contractive complex function with real fixed point
        Inputs = [2.0 + 0.5im]

        result = fixed_point(func, Inputs, Simple(), opts_100)
        @test result.convergence < 1e-8  # Should converge

        if !ismissing(result.fixed_point)
            # Fixed point should be 0.1 / (1 - 0.8) = 0.5
            expected = 0.5 + 0.0im
            @test abs(result.fixed_point[1] - expected) < 1e-8
        end
    end

    @testset "Complex Trigonometric Functions" begin
        # Test with complex trigonometric functions
        func(x) = [0.5 * cos(x[1]) + 0.1im]
        Inputs = [0.5 + 0.1im]

        result = fixed_point(func, Inputs, Anderson(), opts_100)
        @test result.convergence < 1e-8  # Slightly relaxed tolerance for trig functions

        # Verify it's actually a fixed point
        if !ismissing(result.fixed_point)
            residual = abs(func(result.fixed_point)[1] - result.fixed_point[1])
            @test residual < 1e-8
        end
    end

    @testset "Complex Contractive Mapping" begin
        # Test a contractive mapping in the complex plane
        func(x) = [0.7 * x[1] + 0.2 + 0.3im]
        Inputs = [1.0 + 1.0im]

        result = fixed_point(func, Inputs, Simple())
        @test result.convergence < 1e-10

        # The fixed point should be (0.2 + 0.3im) / (1 - 0.7) = (0.2 + 0.3im) / 0.3
        expected_fp = (0.2 + 0.3im) / 0.3
        if !ismissing(result.fixed_point)
            @test abs(result.fixed_point[1] - expected_fp) < 1e-9  # Slightly relaxed tolerance
        end
    end

    @testset "Complex Error Handling" begin
        # Test error handling with complex numbers

        # Function that might produce NaN in complex domain
        problematic_func(x) = [sqrt(x[1])]  # sqrt of negative complex numbers

        # This should handle complex sqrt properly
        result = fixed_point(problematic_func, [-1.0 + 0.0im], Simple(), opts_10)
        # Should not crash - either converge or reach max iterations
        @test result.termination_condition in [
            :ReachedConvergenceThreshold, :ReachedMaxIter, :InvalidInputOrOutputOfIteration
        ]

        # Test with function that returns Inf
        inf_func(x) = [Inf + 0.0im]
        opts_5_quiet = FixedPointOptions(max_iterations=5, quiet_errors=true)
        result_inf = fixed_point(inf_func, [1.0 + 0.0im], Simple(), opts_5_quiet)
        @test result_inf.termination_condition == :InvalidInputOrOutputOfIteration
    end

    @testset "Complex Convergence Metrics" begin
        # Test that convergence metrics work properly with complex numbers
        func(x) = [0.9 * x[1] + 0.05im]
        Inputs = [1.0 + 0.5im]

        result = fixed_point(func, Inputs, Anderson())

        # Convergence should be real-valued even for complex inputs
        @test isa(result.convergence, Real)
        @test result.convergence >= 0  # Convergence metrics should be non-negative

        # ConvergenceVector should contain real values
        if !ismissing(result.convergence_vector)
            @test all(isa.(result.convergence_vector, Real))
            @test all(result.convergence_vector .>= 0)
        end
    end

    @testset "Mixed Real-Complex Operations" begin
        # Test operations that mix real and complex numbers
        func(x) = [ComplexF64(0.5) * x[1] + 0.1]  # Real coefficient, complex result

        # Start with real input
        result1 = fixed_point(func, [1.0], Anderson())
        @test result1.convergence < 1e-10

        # Start with complex input
        result2 = fixed_point(func, [1.0 + 0.0im], Anderson())
        @test result2.convergence < 1e-10

        # Results should be equivalent
        if !ismissing(result1.fixed_point) && !ismissing(result2.fixed_point)
            @test abs(result1.fixed_point[1] - result2.fixed_point[1]) < 1e-10
        end
    end

    @testset "Complex Anderson Acceleration" begin
        # Specific test for Anderson acceleration with complex numbers
        # This tests the pinv() fallback for complex matrices
        func(x) = [0.8 * x[1] + 0.1im, 0.7 * x[2] - 0.05im]
        Inputs = [1.0 + 0.2im, 0.5 - 0.1im]

        result = fixed_point(func, Inputs, Anderson(maxM=3), opts_silent)
        @test result.convergence < 1e-10

        # Test that Anderson method converges faster than simple iteration
        result_simple = fixed_point(func, Inputs, Simple(), opts_100)

        # Anderson should typically use fewer iterations
        if !ismissing(result.fixed_point) && !ismissing(result_simple.fixed_point)
            @test result.iterations <= result_simple.iterations
        end
    end
end

end
