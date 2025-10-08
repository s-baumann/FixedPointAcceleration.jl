
using Test
using FixedPointAccelerationNext
using FixedPointAcceleration

func(x) = cos.(x)
@testset "Simple Scalar Function" begin
    # Core already loaded in runtests; remove duplicate using to avoid warnings

    x0 = [1.1]
    fp_simple = solve(func, x0; method=Simple())
    @test fp_simple.residual_norm < 1e-10
    fp_anderson = solve(func, x0; method=Anderson())
    @test fp_anderson.residual_norm < 1e-10
    fp_aitken = solve(func, x0; method=Aitken())
    @test fp_aitken.residual_norm < 1e-10
    # fp_newton = solve(func, x0; method= Newton())
    # @test fp_newton.residual_norm < 1e-10
    fp_VEA = solve(func, x0; method=VEA())
    @test fp_VEA.residual_norm < 1e-10
    fp_SEA = solve(func, x0; method=SEA())
    @test fp_SEA.residual_norm < 1e-10
    fp_MPE = solve(func, x0; method=MPE())
    @test fp_MPE.residual_norm < 1e-10
    fp_RRE = solve(func, x0; method=RRE())
    @test fp_RRE.residual_norm < 1e-10
end

@testset "new vs old parity" begin
    # Wrap scalar function as vector-valued 1D function for new core
    f_vec(x) = (y=similar(x); y[1]=cos(x[1]); y)
    cfg = FixedPointConfig(;
        threshold=1e-10, max_iters=1_000, relaxation=1.0, relaxation_reference=:Simple
    )
    x0_vec = [1.1]

    legacy_res = fixed_point(func, copy(x0_vec); Algorithm=:Simple)
    legacy_fp_val = legacy_res.FixedPoint_[1]
    # Allow a small numerical slack (methods may reorder history / rounding)
    val_atol = 1e-8
    res_atol = 1e-10

    for (label, method_old, method_new) in [
        ("Simple_parity", :Simple, Simple()),
        ("Anderson_parity", :Anderson, Anderson()),
        ("Aitken_parity", :Aitken, Aitken()),
        ("MPE_parity", :MPE, MPE()),
        ("RRE_parity", :RRE, RRE()),
        ("VEA_parity", :VEA, VEA()),
        ("SEA_parity", :SEA, SEA()),
    ]
        @testset "$label" begin
            sol = solve(f_vec, x0_vec; method=method_new, cfg)
            @test sol.status == :Converged
            @test sol.residual_norm < res_atol
            @test abs(sol.fixed_point[1] - legacy_fp_val) < val_atol
        end
    end

    # In-place variant (vector form)
    f_inplace!(out, x) = (out[1] = cos(x[1]))
    for (label, meth) in [("Simple_inplace", Simple()), ("Anderson_inplace", Anderson())]
        @testset "$label" begin
            sol_inplace = solve!(f_inplace!, copy(x0_vec); method=meth, cfg=cfg)
            @test sol_inplace.status == :Converged
            @test sol_inplace.residual_norm < res_atol
            @test abs(sol_inplace.fixed_point[1] - legacy_fp_val) < val_atol
        end
    end
end
