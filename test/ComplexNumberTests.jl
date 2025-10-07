using Test
@testset "Complex Number Functions" begin
    using FixedPointAcceleration

    @testset "Simple Complex Scalar" begin
        # Test a simple complex function with a known fixed point
        func(x) = [0.5 * x[1] + 0.25im]
        Inputs = [1.0 + 0.5im]
        
        # Test with different algorithms
        fp_simple = fixed_point(func, Inputs; Algorithm=:Simple)
        @test fp_simple.Convergence_ < 1e-10
        @test abs(fp_simple.FixedPoint_[1] - (0.0 + 0.5im)) < 1e-10
        
        fp_anderson = fixed_point(func, Inputs; Algorithm=:Anderson)
        @test fp_anderson.Convergence_ < 1e-10
        @test abs(fp_anderson.FixedPoint_[1] - (0.0 + 0.5im)) < 1e-10
        
        fp_aitken = fixed_point(func, Inputs; Algorithm=:Aitken)
        @test fp_aitken.Convergence_ < 1e-10
        @test abs(fp_aitken.FixedPoint_[1] - (0.0 + 0.5im)) < 1e-10
        
        fp_newton = fixed_point(func, Inputs; Algorithm=:Newton)
        @test fp_newton.Convergence_ < 1e-10
        @test abs(fp_newton.FixedPoint_[1] - (0.0 + 0.5im)) < 1e-10
        
        fp_VEA = fixed_point(func, Inputs; Algorithm=:VEA)
        @test fp_VEA.Convergence_ < 1e-10
        @test abs(fp_VEA.FixedPoint_[1] - (0.0 + 0.5im)) < 1e-10
        
        fp_SEA = fixed_point(func, Inputs; Algorithm=:SEA)
        @test fp_SEA.Convergence_ < 1e-10
        @test abs(fp_SEA.FixedPoint_[1] - (0.0 + 0.5im)) < 1e-10
        
        fp_MPE = fixed_point(func, Inputs; Algorithm=:MPE)
        @test fp_MPE.Convergence_ < 1e-10
        @test abs(fp_MPE.FixedPoint_[1] - (0.0 + 0.5im)) < 1e-10
        
        fp_RRE = fixed_point(func, Inputs; Algorithm=:RRE)
        @test fp_RRE.Convergence_ < 1e-10
        @test abs(fp_RRE.FixedPoint_[1] - (0.0 + 0.5im)) < 1e-10
    end

    @testset "Complex Vector Functions" begin
        # Test a 2D complex function
        func(x) = [0.5 * x[1] + 0.1im, 0.3 * x[2] - 0.2im]
        Inputs = [1.0 + 0.2im, 0.5 - 0.3im]
        
        fp_simple = fixed_point(func, Inputs; Algorithm=:Simple)
        @test fp_simple.Convergence_ < 1e-10
        
        fp_anderson = fixed_point(func, Inputs; Algorithm=:Anderson)
        @test fp_anderson.Convergence_ < 1e-10
        
        # Test that the fixed point satisfies f(x) = x
        if !ismissing(fp_anderson.FixedPoint_)
            residual = func(fp_anderson.FixedPoint_) .- fp_anderson.FixedPoint_
            @test maximum(abs.(residual)) < 1e-10
        end
    end

    @testset "Complex Scalar Input" begin
        # Test with single complex number input (not vector)
        func(x) = [0.8 * x[1] + 0.1im]
        result = fixed_point(func, 1.0 + 0.5im; Algorithm=:Anderson)
        @test result.Convergence_ < 1e-10
        @test !ismissing(result.FixedPoint_)
    end

    @testset "Complex with Real Output" begin
        # Test a function that takes complex input but produces complex output (but real part dominates)
        func(x) = [0.8 * x[1] + 0.1]  # Contractive complex function with real fixed point
        Inputs = [2.0 + 0.5im]
        
        result = fixed_point(func, Inputs; Algorithm=:Simple, MaxIter=100)
        @test result.Convergence_ < 1e-8  # Should converge
        
        if !ismissing(result.FixedPoint_)
            # Fixed point should be 0.1 / (1 - 0.8) = 0.5
            expected = 0.5 + 0.0im
            @test abs(result.FixedPoint_[1] - expected) < 1e-8
        end
    end

    @testset "Complex Trigonometric Functions" begin
        # Test with complex trigonometric functions
        func(x) = [0.5 * cos(x[1]) + 0.1im]
        Inputs = [0.5 + 0.1im]
        
        result = fixed_point(func, Inputs; Algorithm=:Anderson, MaxIter=100)
        @test result.Convergence_ < 1e-8  # Slightly relaxed tolerance for trig functions
        
        # Verify it's actually a fixed point
        if !ismissing(result.FixedPoint_)
            residual = abs(func(result.FixedPoint_)[1] - result.FixedPoint_[1])
            @test residual < 1e-8
        end
    end

    @testset "Complex Contractive Mapping" begin
        # Test a contractive mapping in the complex plane
        func(x) = [0.7 * x[1] + 0.2 + 0.3im]
        Inputs = [1.0 + 1.0im]
        
        result = fixed_point(func, Inputs; Algorithm=:Simple)
        @test result.Convergence_ < 1e-10
        
        # The fixed point should be (0.2 + 0.3im) / (1 - 0.7) = (0.2 + 0.3im) / 0.3
        expected_fp = (0.2 + 0.3im) / 0.3
        if !ismissing(result.FixedPoint_)
            @test abs(result.FixedPoint_[1] - expected_fp) < 1e-9  # Slightly relaxed tolerance
        end
    end

    @testset "Complex Error Handling" begin
        # Test error handling with complex numbers
        
        # Function that might produce NaN in complex domain
        problematic_func(x) = [sqrt(x[1])]  # sqrt of negative complex numbers
        
        # This should handle complex sqrt properly
        result = fixed_point(problematic_func, [-1.0 + 0.0im]; Algorithm=:Simple, MaxIter=10)
        # Should not crash - either converge or reach max iterations
        @test result.TerminationCondition_ in [:ReachedConvergenceThreshold, :ReachedMaxIter, :InvalidInputOrOutputOfIteration]
        
        # Test with function that returns complex infinity
        inf_func(x) = [x[1] * (1.0 + 0.0im) / 0.0]
        result_inf = fixed_point(inf_func, [1.0 + 0.0im]; Algorithm=:Simple, MaxIter=5, quiet_errors=true)
        @test result_inf.TerminationCondition_ == :InvalidInputOrOutputOfIteration
    end

    @testset "Complex Convergence Metrics" begin
        # Test that convergence metrics work properly with complex numbers
        func(x) = [0.9 * x[1] + 0.05im]
        Inputs = [1.0 + 0.5im]
        
        result = fixed_point(func, Inputs; Algorithm=:Anderson)
        
        # Convergence should be real-valued even for complex inputs
        @test isa(result.Convergence_, Real)
        @test result.Convergence_ >= 0  # Convergence metrics should be non-negative
        
        # ConvergenceVector should contain real values
        if !ismissing(result.ConvergenceVector_)
            @test all(isa.(result.ConvergenceVector_, Real))
            @test all(result.ConvergenceVector_ .>= 0)
        end
    end

    @testset "Mixed Real-Complex Operations" begin
        # Test operations that mix real and complex numbers
        func(x) = [ComplexF64(0.5) * x[1] + 0.1]  # Real coefficient, complex result
        
        # Start with real input
        result1 = fixed_point(func, [1.0]; Algorithm=:Anderson)
        @test result1.Convergence_ < 1e-10
        
        # Start with complex input  
        result2 = fixed_point(func, [1.0 + 0.0im]; Algorithm=:Anderson)
        @test result2.Convergence_ < 1e-10
        
        # Results should be equivalent
        if !ismissing(result1.FixedPoint_) && !ismissing(result2.FixedPoint_)
            @test abs(result1.FixedPoint_[1] - result2.FixedPoint_[1]) < 1e-10
        end
    end

    @testset "Complex Anderson Acceleration" begin
        # Specific test for Anderson acceleration with complex numbers
        # This tests the pinv() fallback for complex matrices
        func(x) = [0.8 * x[1] + 0.1im, 0.7 * x[2] - 0.05im]
        Inputs = [1.0 + 0.2im, 0.5 - 0.1im]
        
        result = fixed_point(func, Inputs; Algorithm=:Anderson, MaxM=3, PrintReports=false)
        @test result.Convergence_ < 1e-10
        
        # Test that Anderson method converges faster than simple iteration
        result_simple = fixed_point(func, Inputs; Algorithm=:Simple, MaxIter=100)
        
        # Anderson should typically use fewer iterations
        if !ismissing(result.FixedPoint_) && !ismissing(result_simple.FixedPoint_)
            @test result.Iterations_ <= result_simple.Iterations_
        end
    end
end