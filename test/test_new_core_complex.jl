
using Test, FixedPointAccelerationNext
using FixedPointAccelerationNext:
    solve, solve!, FixedPointConfig, Simple, Anderson, Aitken, MPE, RRE, VEA, SEA
using LinearAlgebra: norm

@testset "complex scalar contraction" begin
    α = 0.5 + 0.2im
    β = 0.3 - 0.1im
    f(x) = α .* x .+ β
    x0 = [1.0 + 0.5im]
    expected = β / (1 - α)
    cfg = FixedPointConfig(; threshold=1e-11, max_iters=120)
    sol_simple = solve(f, x0; method=Simple(), cfg=cfg)
    @test sol_simple.status == :Converged
    @test sol_simple.fixed_point[1] ≈ expected atol=1e-9

    sol_anderson = solve(f, x0; method=Anderson(m=4), cfg=cfg)
    @test sol_anderson.fixed_point[1] ≈ expected atol=1e-9
    @test sol_anderson.iterations <= sol_simple.iterations

    sol_aitken = solve(f, x0; method=Aitken(), cfg=cfg)
    @test sol_aitken.fixed_point[1] ≈ expected atol=1e-9
    @test sol_aitken.iterations <= sol_simple.iterations

    sol_mpe = solve(f, x0; method=MPE(period=3), cfg=cfg)
    @test sol_mpe.fixed_point[1] ≈ expected atol=1e-9
    @test sol_mpe.iterations <= sol_simple.iterations

    sol_rre = solve(f, x0; method=RRE(period=3), cfg=cfg)
    @test sol_rre.fixed_point[1] ≈ expected atol=1e-9
    @test sol_rre.iterations <= sol_simple.iterations

    sol_vea = solve(f, x0; method=VEA(period=5), cfg=cfg)
    @test sol_vea.fixed_point[1] ≈ expected atol=1e-9
    @test sol_vea.iterations <= sol_simple.iterations

    sol_sea = solve(f, x0; method=SEA(period=6), cfg=cfg)
    @test sol_sea.fixed_point[1] ≈ expected atol=1e-9
    @test sol_sea.iterations <= sol_simple.iterations
end

@testset "complex vector contraction" begin
    α1 = 0.4 + 0.15im
    β1 = 0.2 - 0.05im
    α2 = 0.3 - 0.25im
    β2 = -0.1 + 0.2im
    f(x) = [α1 * x[1] + β1, α2 * x[2] + β2]
    x0 = [1.0 + 0.3im, -0.5 + 0.7im]
    expected = [(β1) / (1 - α1), (β2) / (1 - α2)]
    cfg = FixedPointConfig(; threshold=1e-11, max_iters=150)
    sol_simple = solve(f, x0; method=Simple(), cfg=cfg)
    @test sol_simple.status == :Converged
    @test sol_simple.fixed_point ≈ expected atol=1e-9
    sol_anderson = solve(f, x0; method=Anderson(m=5), cfg=cfg)
    @test sol_anderson.fixed_point ≈ expected atol=1e-9
    @test sol_anderson.iterations <= sol_simple.iterations
end

@testset "complex residual norm real" begin
    f(x) = 0.6 .* x .+ (0.1 + 0.05im)
    x0 = [2.0 + 1.0im]
    cfg = FixedPointConfig(; threshold=1e-12)
    sol = solve(f, x0; method=Simple(), cfg=cfg)
    @test sol.status == :Converged
    @test isa(sol.residual_norm, Real)
    @test sol.residual_norm >= 0
end

@testset "in-place complex solve!" begin
    α = 0.45 - 0.1im
    β = -0.2 + 0.3im
    function f!(out, x)
        @. out = α * x + β
    end
    x = [1.2 + 0.8im]
    cfg = FixedPointConfig(; threshold=1e-12, max_iters=120)
    sol = solve!(f!, x; method=Anderson(m=3), cfg=cfg)
    expected = β / (1 - α)
    @test sol.status == :Converged
    @test sol.fixed_point[1] ≈ expected atol=1e-9
end

@testset "complex high-dim polynomial & epsilon" begin
    # Construct element-wise complex contraction mapping of dimension 5
    α = ComplexF64[
        0.55 + 0.10im, 0.40 - 0.20im, 0.30 + 0.25im, 0.65 - 0.05im, 0.50 + 0.30im
    ]
    β = ComplexF64[-0.2 + 0.3im, 0.15 - 0.10im, 0.05 + 0.05im, -0.1 + 0.2im, 0.25 - 0.15im]
    f(x) = α .* x .+ β
    x0 = [0.8 - 0.4im, -0.3 + 0.9im, 1.1 - 0.2im, -0.5 + 0.7im, 0.6 + 0.3im]
    expected = β ./ (1 .- α)
    # Looser threshold and higher iteration cap: mixed complex coefficients can
    # reduce effectiveness of generic accelerators (esp. without specialized LS solves).
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=400)

    base = solve(f, x0; method=Simple(), cfg=cfg)
    @test base.status == :Converged
    @test base.fixed_point ≈ expected atol=1e-8

    # Helper to evaluate closeness via norm and ensure not worse than baseline residual by >10%
    function _assert_close(name, sol)
        err = norm(sol.fixed_point .- expected)
        @test err < 1e-6
        @test sol.residual_norm <= 1.10 * base.residual_norm
    end

    mpe = solve(f, x0; method=MPE(period=3), cfg=cfg);
    _assert_close("MPE", mpe)
    rre = solve(f, x0; method=RRE(period=3), cfg=cfg);
    _assert_close("RRE", rre)
    vea = solve(f, x0; method=VEA(period=5), cfg=cfg);
    _assert_close("VEA", vea)
    sea = solve(f, x0; method=SEA(period=6), cfg=cfg);
    _assert_close("SEA", sea)
end
