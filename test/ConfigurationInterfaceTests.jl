using FixedPointAcceleration
using Test

@testset "Configuration Interface Tests" begin
    # Test function for simple vector case
    g(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]

    @testset "Configuration Struct Creation" begin
        # Test default construction
        opts_default = default_options()
        @test opts_default isa FixedPointOptions
        @test opts_default.convergence.threshold == 1e-10
        @test opts_default.stability.dampening == 1.0
        @test opts_default.reporting.print_reports == false

        # Test preset configurations
        opts_robust = robust_options()
        @test opts_robust.stability.dampening == 0.7
        @test opts_robust.stability.replace_invalids == :ReplaceElements
        @test opts_robust.stability.quiet_errors == true

        opts_fast = fast_options()
        @test opts_fast.convergence.threshold == 1e-12
        @test opts_fast.convergence.max_iterations == 500
        @test opts_fast.stability.dampening == 1.0

        opts_debug = debug_options()
        @test opts_debug.reporting.print_reports == true
        @test opts_debug.reporting.reporting_sig_figs == 12
        @test opts_debug.stability.replace_invalids == :ReplaceElements

        # Test custom configuration
        custom_opts = FixedPointOptions(
            threshold = 1e-8,
            max_iterations = 50,
            dampening = 0.8,
            print_reports = false
        )
        @test custom_opts.convergence.threshold == 1e-8
        @test custom_opts.convergence.max_iterations == 50
        @test custom_opts.stability.dampening == 0.8
        @test custom_opts.reporting.print_reports == false
    end

    @testset "New Fixed Point Interface" begin
        # Test new interface with default options
        result1 = fixed_point(g, [0.3, 0.9], Anderson(), default_options())
        @test result1 isa FixedPointResults
        @test !ismissing(result1.FixedPoint_)
        @test result1.TerminationCondition_ == :ReachedConvergenceThreshold

        # Test with preset options
        result2 = fixed_point(g, [0.3, 0.9], Anderson(), robust_options())
        @test result2 isa FixedPointResults
        @test !ismissing(result2.FixedPoint_)

        # Test with custom options
        custom_opts = FixedPointOptions(threshold=1e-8, max_iterations=30)
        result3 = fixed_point(g, [0.3, 0.9], Anderson(), custom_opts)
        @test result3 isa FixedPointResults
        @test result3.Iterations_ <= 30  # Should respect max_iterations

        # Test continuing from previous results
        result4 = fixed_point(g, result1, Simple(), robust_options())
        @test result4 isa FixedPointResults
        @test result4.Iterations_ >= result1.Iterations_  # Should have more iterations
    end

    @testset "Backward Compatibility" begin
        # Test that old interface still works
        opts_old_equiv = FixedPointOptions(max_iterations=20, print_reports=false)
        result_old = fixed_point(g, [0.3, 0.9], Anderson(), opts_old_equiv)
        @test result_old isa FixedPointResults
        @test !ismissing(result_old.FixedPoint_)
        @test result_old.Iterations_ <= 20

        # Results should be equivalent between old and new interfaces
        result_new = fixed_point(g, [0.3, 0.9], Anderson(),
                                FixedPointOptions(max_iterations=20, print_reports=false))

        # Should get same convergence (within numerical tolerance)
        @test isapprox(result_old.FixedPoint_, result_new.FixedPoint_, rtol=1e-10)
    end
end
