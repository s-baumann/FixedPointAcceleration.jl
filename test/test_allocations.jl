module test_allocations

using Test
using FixedPointAcceleration: solve, FixedPointConfig, Anderson, Simple, MPE, RRE, VEA, SEA
using LinearAlgebra: norm

# NOTE: These tests are intentionally coarse. They serve as regression guards
# to catch large accidental allocation increases in the new core algorithms.
# Absolute numbers may vary slightly across Julia versions / BLAS libs, so
# thresholds are padded. If these start failing spuriously, adjust upward
# conservatively, but investigate first.

# Helper to warm-up JIT before measuring allocations
function warmup_and_measure(f)
    f() # warmup
    GC.gc()
    return @allocated f()
end

@testset "vea allocation amortization" begin
    f(x) = 0.57 .* x .+ 0.9
    x0 = rand(10)
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=45, history_window=14)
    method = VEA(period=5)
    solve(f, x0; method=method, cfg=cfg) # warmup (JIT + allocations)
    GC.gc()
    alloc_second = @allocated solve(f, x0; method=method, cfg=cfg)
    # Coarse ceiling; epsilon currently uses fallback (allocating), keep generous.
    @test alloc_second < 650_000
end

@testset "sea allocation amortization" begin
    f(x) = 0.61 .* x .+ 0.7
    x0 = rand(10)
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=45, history_window=14)
    method = SEA(period=6)
    solve(f, x0; method=method, cfg=cfg)
    GC.gc()
    alloc_second = @allocated solve(f, x0; method=method, cfg=cfg)
    @test alloc_second < 550_000
end

@testset "mpe allocation amortization" begin
    f(x) = 0.55 .* x .+ 0.8
    x0 = rand(10)
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=40, history_window=10)
    method = MPE(period=3)
    # warmup
    solve(f, x0; method=method, cfg=cfg)
    GC.gc()
    alloc_second = @allocated solve(f, x0; method=method, cfg=cfg)
    @test alloc_second < 450_000
end

@testset "rre allocation amortization" begin
    f(x) = 0.52 .* x .+ 1.2
    x0 = rand(10)
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=50, history_window=12)
    method = RRE(period=4)
    solve(f, x0; method=method, cfg=cfg)
    GC.gc()
    alloc_second = @allocated solve(f, x0; method=method, cfg=cfg)
    @test alloc_second < 500_000
end

@testset "allocation simple vs anderson (small problem)" begin
    f(x) = 0.5 .* x .+ 1.0
    x0 = fill(3.0, 5)
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=20)

    alloc_simple = warmup_and_measure() do
        solve(f, x0; method=Simple(), cfg=cfg)
    end

    alloc_anderson = warmup_and_measure() do
        solve(f, x0; method=Anderson(m=4), cfg=cfg)
    end

    # Anderson will allocate more due to its linear solves, but both should be
    # comfortably below these coarse upper bounds.
    @test alloc_simple < 200_000
    @test alloc_anderson < 600_000
end

@testset "anderson repeated call amortization" begin
    f(x) = 0.3 .* x .+ 2.5
    x0 = rand(8)
    cfg = FixedPointConfig(; threshold=1e-9, max_iters=25, history_window=5)
    method = Anderson(m=5)

    # First run (JIT + allocations to populate workspace)
    solve(f, x0; method=method, cfg=cfg)
    GC.gc()
    alloc_second = @allocated solve(f, x0; method=method, cfg=cfg)

    # After warm-up, expect materially reduced allocations. Set a soft ceiling.
    @test alloc_second < 400_000
end

@testset "simple method stable low allocations on repeat" begin
    f(x) = 0.4 .* x .+ 1.0
    x0 = rand(6)
    cfg = FixedPointConfig(; threshold=1e-9, max_iters=30)

    solve(f, x0; method=Simple(), cfg=cfg) # warmup
    GC.gc()
    alloc_second = @allocated solve(f, x0; method=Simple(), cfg=cfg)
    @test alloc_second < 120_000
end

end # module test_allocations
