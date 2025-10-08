module test_simple_scalar

using Test
using FixedPointAcceleration.OldImplementation:
    fixed_point,
    Simple as LegacySimple,
    Anderson as LegacyAnderson,
    Aitken as LegacyAitken,
    Newton as LegacyNewton,
    MPE as LegacyMPE,
    RRE as LegacyRRE,
    VEA as LegacyVEA,
    SEA as LegacySEA,
    FixedPointOptions
@testset "Simple Scalar Function" begin
    # Core already loaded in runtests; remove duplicate using to avoid warnings

    func(x) = cos.(x)
    Inputs = 1.1
    fp_simple = fixed_point(func, Inputs, LegacySimple())
    @test fp_simple.convergence < 1e-10
    fp_anderson = fixed_point(func, Inputs, LegacyAnderson())
    @test fp_anderson.convergence < 1e-10
    fp_aitken = fixed_point(func, Inputs, LegacyAitken())
    @test fp_aitken.convergence < 1e-10
    fp_newton = fixed_point(func, Inputs, LegacyNewton())
    @test fp_newton.convergence < 1e-10
    fp_VEA = fixed_point(func, Inputs, LegacyVEA())
    @test fp_VEA.convergence < 1e-10
    fp_SEA = fixed_point(func, Inputs, LegacySEA())
    @test fp_SEA.convergence < 1e-10
    fp_MPE = fixed_point(func, Inputs, LegacyMPE())
    @test fp_MPE.convergence < 1e-10
    fp_RRE = fixed_point(func, Inputs, LegacyRRE())
    @test fp_RRE.convergence < 1e-10

    # Testing input/output matrix functionality
    opts_20 = FixedPointOptions(max_iterations=20)
    a = fixed_point(func, fp_RRE.inputs, fp_RRE.outputs, LegacyRRE(), opts_20)
    @test a.convergence < 1e-8 # If the above line throws then that is a test failure.
    # Testing different input/output shapes
    b = fixed_point(
        func, fp_RRE.inputs[:, 1:1], fp_RRE.outputs[:, 1:1], LegacyRRE(), opts_20
    )
    @test b.convergence < 1e-8 # If the above line throws then that is a test failure.

    # === New core (solve) parity tests ===
    using FixedPointAcceleration:
        FixedPointConfig, solve, solve!, Simple, Anderson, Aitken, MPE, RRE, VEA, SEA
    # Wrap scalar function as vector-valued 1D function for new core
    f_vec(x) = (y=similar(x); y[1]=cos(x[1]); y)
    cfg = FixedPointConfig(;
        threshold=1e-10, max_iters=1_000, relaxation=1.0, relaxation_reference=:Simple
    )
    x0_vec = [Inputs]

    legacy_fp_val = fp_simple.fixed_point[1]
    # Allow a small numerical slack (methods may reorder history / rounding)
    val_atol = 1e-8
    res_atol = 1e-10

    for (label, meth_old, meth_new) in [
        ("Simple_parity", LegacySimple(), Simple()),
        ("Anderson_parity", LegacyAnderson(), Anderson()),
        ("Aitken_parity", LegacyAitken(), Aitken()),
        ("MPE_parity", LegacyMPE(), MPE()),
        ("RRE_parity", LegacyRRE(), RRE()),
        ("VEA_parity", LegacyVEA(), VEA()),
        ("SEA_parity", LegacySEA(), SEA()),
    ]
        @testset "$label" begin
            sol = solve(f_vec, x0_vec; method=meth_new, cfg=cfg)
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
end
