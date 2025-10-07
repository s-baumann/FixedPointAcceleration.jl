module test_new_core

using Test
using FixedPointAcceleration:
    solve, solve!, FixedPointConfig, Anderson, Simple, Aitken, MPE, RRE, VEA, SEA
using LinearAlgebra: norm

@testset "simple linear contraction" begin
    f(x) = 0.5 .* x .+ 1.0
    x0 = [10.0, -4.0]
    sol = solve(f, x0; method=Simple(), cfg=FixedPointConfig(; threshold=1e-8))
    @test sol.status == :Converged
    @test norm(sol.fixed_point .- (2.0 .* ones(2))) < 1e-6
end

@testset "anderson matches simple fixed point" begin
    f(x) = 0.3 .* x .+ 3.0
    x0 = [0.0, 0.0]
    sol_anderson = solve(
        f, x0; method=Anderson(m=4), cfg=FixedPointConfig(; threshold=1e-9)
    )
    sol_simple = solve(f, x0; method=Simple(), cfg=FixedPointConfig(; threshold=1e-9))
    @test sol_anderson.fixed_point ≈ sol_simple.fixed_point atol=1e-9
end

@testset "aitken acceleration" begin
    f(x) = 0.5 .* x .+ 1.0
    x0 = [10.0]
    sol_simple = solve(f, x0; method=Simple(), cfg=FixedPointConfig(; threshold=1e-12))
    sol_aitken = solve(f, x0; method=Aitken(), cfg=FixedPointConfig(; threshold=1e-12))
    @test sol_simple.fixed_point ≈ sol_aitken.fixed_point atol=1e-10
    @test sol_aitken.iterations <= sol_simple.iterations
end

@testset "divergence detection" begin
    f(x) = 2.0 .* x
    x0 = [1.0]
    cfg = FixedPointConfig(; max_iters=50, divergence_factor=10.0, threshold=1e-30)
    sol = solve(f, x0; method=Simple(), cfg=cfg)
    @test sol.status == :Diverged
end

@testset "history window enforcement" begin
    f(x) = 0.5 .* x
    x0 = [4.0]
    cfg = FixedPointConfig(; threshold=1e-6, history_window=3)
    sol = solve(f, x0; method=Simple(), cfg=cfg)
    @test length(sol.history) <= 3
end

@testset "in-place solve!" begin
    function f!(out, x)
        @. out = 0.25 * x + 3
    end
    x = [0.0]
    sol = solve!(f!, x; method=Simple(), cfg=FixedPointConfig(; threshold=1e-12))
    @test sol.status == :Converged
    @test sol.fixed_point[1] ≈ 4.0 atol=1e-8
end

@testset "mpe acceleration periodic" begin
    f(x) = 0.6 .* x .+ 1.0
    x0 = [5.0]
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=200)
    sol_simple = solve(f, x0; method=Simple(), cfg=cfg)
    sol_mpe = solve(f, x0; method=MPE(period=3), cfg=cfg)
    @test sol_mpe.fixed_point ≈ sol_simple.fixed_point atol=1e-8
    @test sol_mpe.iterations <= sol_simple.iterations
end

@testset "rre acceleration periodic" begin
    f(x) = 0.55 .* x .+ 2.0
    x0 = [0.0]
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=200)
    sol_simple = solve(f, x0; method=Simple(), cfg=cfg)
    sol_rre = solve(f, x0; method=RRE(period=4), cfg=cfg)
    @test sol_rre.fixed_point ≈ sol_simple.fixed_point atol=1e-8
    @test sol_rre.iterations <= sol_simple.iterations
end

@testset "vea acceleration periodic" begin
    f(x) = 0.62 .* x .+ 0.7
    x0 = [3.0]
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=220)
    sol_simple = solve(f, x0; method=Simple(), cfg=cfg)
    sol_vea = solve(f, x0; method=VEA(period=5), cfg=cfg)
    @test sol_vea.fixed_point ≈ sol_simple.fixed_point atol=1e-8
    @test sol_vea.iterations <= sol_simple.iterations
end

@testset "sea acceleration periodic" begin
    f(x) = 0.58 .* x .+ 1.3
    x0 = [0.5]
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=220)
    sol_simple = solve(f, x0; method=Simple(), cfg=cfg)
    sol_sea = solve(f, x0; method=SEA(period=6), cfg=cfg)
    @test sol_sea.fixed_point ≈ sol_simple.fixed_point atol=1e-8
    @test sol_sea.iterations <= sol_simple.iterations
end

end # module test_new_core
