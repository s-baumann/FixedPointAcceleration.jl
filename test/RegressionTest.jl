
# Regression test for FixedPointAcceleration.jl
# This file checks that numerical results remain stable across code changes.
# If an algorithm's convergence behaviour changes (different iteration count,
# different fixed-point value, or different convergence metric), these tests
# will fail, alerting us to unintended numerical differences.

using FixedPointAcceleration, Test

# ──────────────────────────────────────────────────────────────────────────────
# Test functions
# ──────────────────────────────────────────────────────────────────────────────
const DOTTIE = 0.7390851332151607  # cos(x) = x

func_scalar(x) = cos.(x)
func_vector(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
func_complex_scalar(z) = (z .+ (1 + 2im)) ./ 2
func_complex_vector(v) = (v .+ ComplexF64[1 + 2im, -3 + 1im]) ./ 2
func_sqrt(x) = sqrt.(x)

# ──────────────────────────────────────────────────────────────────────────────
# Reference values: (algorithm, iterations, fixed_point, convergence)
# Captured from v1.1.0 on Julia 1.x
# ──────────────────────────────────────────────────────────────────────────────

@testset "Regression Tests" begin

    # ── Scalar: cos(x), x0 = 1.1 ────────────────────────────────────────────
    @testset "Scalar cos(x)" begin
        scalar_refs = [
            (:Simple,   59, [DOTTIE], 1e-10),
            (:Anderson,  7, [DOTTIE], 1e-10),
            (:Aitken,   10, [DOTTIE], 1e-10),
            (:Newton,    8, [DOTTIE], 1e-10),
            (:VEA,      13, [DOTTIE], 1e-10),
            (:SEA,      13, [DOTTIE], 1e-10),
            (:MPE,      19, [DOTTIE], 1e-10),
            (:RRE,      25, [DOTTIE], 1e-10),
        ]
        for (alg, ref_iters, ref_fp, conv_tol) in scalar_refs
            r = fixed_point(func_scalar, 1.1; Algorithm = alg)
            @testset "$alg" begin
                @test r.TerminationCondition_ == :ReachedConvergenceThreshold
                @test r.Iterations_ == ref_iters
                @test r.FixedPoint_ ≈ ref_fp atol = conv_tol
                @test r.Convergence_ < conv_tol
            end
        end
    end

    # ── Vector: [0.5√|x₁+x₂|, 1.5x₁+0.5x₂], x0 = [0.3, 900.0] ───────────
    @testset "Vector function" begin
        vector_refs = [
            (:Simple,   105, [1.0, 3.0], 1e-9),
            (:Anderson,  12, [1.0, 3.0], 1e-10),
            (:Aitken,   128, [1.0, 3.0], 1e-9),
            (:Newton,   108, [1.0, 3.0], 1e-9),
            (:VEA,       20, [1.0, 3.0], 1e-9),
            (:SEA,       25, [1.0, 3.0], 1e-10),
            (:MPE,       31, [1.0, 3.0], 1e-10),
            (:RRE,       31, [1.0, 3.0], 1e-10),
        ]
        for (alg, ref_iters, ref_fp, conv_tol) in vector_refs
            r = fixed_point(func_vector, [0.3, 900.0]; Algorithm = alg)
            @testset "$alg" begin
                @test r.TerminationCondition_ == :ReachedConvergenceThreshold
                @test r.Iterations_ == ref_iters
                @test r.FixedPoint_ ≈ ref_fp atol = conv_tol
                @test r.Convergence_ < 1e-10
            end
        end
    end

    # ── Complex scalar: (z+(1+2im))/2, z0 = 0+0im ──────────────────────────
    @testset "Complex scalar" begin
        complex_scalar_refs = [
            (:Simple,   35, ComplexF64[1 + 2im], 1e-9),
            (:Anderson,  3, ComplexF64[1 + 2im], 1e-14),
            (:Aitken,    4, ComplexF64[1 + 2im], 1e-14),
            (:Newton,    4, ComplexF64[1 + 2im], 1e-14),
            (:VEA,       7, ComplexF64[1 + 2im], 1e-14),
            (:SEA,       7, ComplexF64[1 + 2im], 1e-14),
            (:MPE,       7, ComplexF64[1 + 2im], 1e-14),
            (:RRE,       7, ComplexF64[1 + 2im], 1e-14),
        ]
        for (alg, ref_iters, ref_fp, conv_tol) in complex_scalar_refs
            r = fixed_point(func_complex_scalar, 0.0 + 0.0im; Algorithm = alg, MaxIter = 200)
            @testset "$alg" begin
                @test r.TerminationCondition_ == :ReachedConvergenceThreshold
                @test r.Iterations_ == ref_iters
                @test r.FixedPoint_ ≈ ref_fp atol = conv_tol
                @test r.Convergence_ < 1e-10
            end
        end
    end

    # ── Complex vector: (v+[1+2im,-3+1im])/2, v0 = [0.25-0.5im, 0.75+0.1im]
    @testset "Complex vector" begin
        complex_vector_refs = [
            (:Simple,   36, ComplexF64[1 + 2im, -3 + 1im], 1e-9),
            (:Anderson,  3, ComplexF64[1 + 2im, -3 + 1im], 1e-14),
            (:Aitken,    4, ComplexF64[1 + 2im, -3 + 1im], 1e-14),
            (:Newton,    4, ComplexF64[1 + 2im, -3 + 1im], 1e-14),
            (:VEA,       7, ComplexF64[1 + 2im, -3 + 1im], 1e-14),
            (:SEA,       7, ComplexF64[1 + 2im, -3 + 1im], 1e-14),
            (:MPE,       7, ComplexF64[1 + 2im, -3 + 1im], 1e-14),
            (:RRE,       7, ComplexF64[1 + 2im, -3 + 1im], 1e-14),
        ]
        for (alg, ref_iters, ref_fp, conv_tol) in complex_vector_refs
            r = fixed_point(func_complex_vector, ComplexF64[0.25 - 0.5im, 0.75 + 0.1im]; Algorithm = alg, MaxIter = 200)
            @testset "$alg" begin
                @test r.TerminationCondition_ == :ReachedConvergenceThreshold
                @test r.Iterations_ == ref_iters
                @test r.FixedPoint_ ≈ ref_fp atol = conv_tol
                @test r.Convergence_ < 1e-10
            end
        end
    end

    # ── Higher-dimensional: sqrt.(x), x0 = [2,3,4,5,7,9,25.6] ──────────────
    @testset "Sqrt vector (7D)" begin
        ref_fp = ones(7)
        sqrt_refs = [
            (:Simple,   :ReachedConvergenceThreshold),
            (:Anderson, :ReachedConvergenceThreshold),
            (:Aitken,   :ReachedConvergenceThreshold),
            (:SEA,      :ReachedConvergenceThreshold),
            (:MPE,      :ReachedConvergenceThreshold),
            (:RRE,      :ReachedConvergenceThreshold),
        ]
        for (alg, expected_term) in sqrt_refs
            r = fixed_point(func_sqrt, [2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 25.6]; Algorithm = alg)
            @testset "$alg" begin
                @test r.TerminationCondition_ == expected_term
                @test r.FixedPoint_ ≈ ref_fp atol = 1e-9
                @test r.Convergence_ < 1e-10
            end
        end
    end

    # ── Dampening variants ───────────────────────────────────────────────────
    @testset "Dampening" begin
        @testset "Anderson Dampening=0.5 (output)" begin
            r = fixed_point(func_scalar, 1.1; Algorithm = :Anderson, Dampening = 0.5)
            @test r.TerminationCondition_ == :ReachedConvergenceThreshold
            @test r.Iterations_ == 35
            @test r.FixedPoint_ ≈ [DOTTIE] atol = 1e-10
        end

        @testset "Anderson Dampening=0.5 (input)" begin
            r = fixed_point(func_scalar, 1.1; Algorithm = :Anderson, Dampening = 0.5, Dampening_With_Input = true)
            @test r.TerminationCondition_ == :ReachedConvergenceThreshold
            @test r.Iterations_ == 21
            @test r.FixedPoint_ ≈ [DOTTIE] atol = 1e-10
        end
    end

    # ── Chaining: pass results from one algorithm into another ───────────────
    @testset "Chaining" begin
        r1 = fixed_point(func_scalar, 1.1; Algorithm = :Simple, MaxIter = 2)
        r2 = fixed_point(func_scalar, r1; Algorithm = :Anderson)
        @test r2.TerminationCondition_ == :ReachedConvergenceThreshold
        @test r2.FixedPoint_ ≈ [DOTTIE] atol = 1e-10
        # Chained result should use more total iterations than Anderson alone
        @test (r1.Iterations_ + r2.Iterations_) > 7
    end

end
