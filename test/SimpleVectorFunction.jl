module test_simple_vector

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
# Core already loaded in runtests; remove duplicate using to avoid warnings
@testset "Simple Vector Function" begin
    func(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
    Inputs = [0.3, 900.0]
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

    # Now trying a function with a typeswitch from Int to Float
    Inputs = [1, 4]
    fp_simple = fixed_point(func, Inputs, LegacySimple())
    @test fp_simple.convergence < 1e-10
    fp_anderson = fixed_point(func, Inputs, LegacyAnderson())
    @test fp_anderson.convergence < 1e-10

    # And finally printing a couple of them to test printing.
    print_opts = FixedPointOptions(print_reports=false)
    fp_anderson2 = fixed_point(func, Inputs, LegacyAnderson(), print_opts)
    @test fp_anderson2.convergence < 1e-10
    fp_simple2 = fixed_point(func, Inputs, LegacySimple(), print_opts)
    @test fp_simple2.convergence < 1e-10

    # Testing function that returns a missing
    f1(x) = [missing, missing]
    fp_simple3 = fixed_point(f1, Inputs, LegacySimple(), print_opts)
    @test fp_simple3.termination_condition == :InvalidInputOrOutputOfIteration
    f2(x) = [1, missing]
    fp_simple4 = fixed_point(f2, Inputs, LegacySimple(), print_opts)
    @test fp_simple4.termination_condition == :InvalidInputOrOutputOfIteration

    # Testing the outputting of a tuple
    f3(x) = (a=[1.0, 1.0], b=:goodProgressInFunction)
    @test_throws ErrorException(
        "This function returned a @NamedTuple{a::Vector{Float64}, b::Symbol}. The Fixedpoint function can only return a vector or a tuple of which the first entry is the vector for which a fixedpoint is sought and the second is a namedtuple (the contents of which are output for the user but are not used in fixed point acceleration).",
    ) fixed_point(f3, Inputs, LegacySimple(), print_opts)
    # And doing side effects properly.
    f4(x) = ([1.0, 1.0], (b=:goodProgressInFunction,))
    fp_simple5 = fixed_point(f4, Inputs, LegacySimple(), print_opts)
    @test fp_simple5.convergence < 1e-10

    # === New core (solve) parity tests (vector) ===
    using FixedPointAcceleration: FixedPointConfig, solve, solve!, Simple, Anderson, Aitken, MPE, RRE, VEA, SEA
    # Configuration for new core
    cfg = FixedPointConfig(;
        threshold=1e-10, max_iters=2_000, relaxation=1.0, relaxation_reference=:Simple
    )
    # Use original float input vector as baseline
    parity_inputs = [0.3, 900.0]
    legacy_simple_vec = fixed_point(func, parity_inputs, LegacySimple()).fixed_point
    val_atol = 1e-8
    res_atol = 1e-10

    for (label, meth_legacy, meth_new) in [
        ("Simple_parity", LegacySimple(), Simple()),
        ("Anderson_parity", LegacyAnderson(), Anderson()),
        ("Aitken_parity", LegacyAitken(), Aitken()),
        ("MPE_parity", LegacyMPE(), MPE()),
        ("RRE_parity", LegacyRRE(), RRE()),
        ("VEA_parity", LegacyVEA(), VEA()),
        ("SEA_parity", LegacySEA(), SEA()),
    ]
        @testset "vector_parity_$label" begin
            sol = solve(func, copy(parity_inputs); method=meth_new, cfg=cfg)
            @test sol.status == :Converged
            @test sol.residual_norm < res_atol
            @test maximum(abs.(sol.fixed_point .- legacy_simple_vec)) < val_atol
        end
    end

    # In-place variants for a subset (Simple, Anderson)
    func_inplace!(out, x) = (out[1]=0.5*sqrt(abs(x[1] + x[2])); out[2]=1.5*x[1] + 0.5*x[2])
    for (label, meth) in [("Simple_inplace", Simple()), ("Anderson_inplace", Anderson())]
        @testset "vector_parity_$label" begin
            sol_inplace = solve!(func_inplace!, copy(parity_inputs); method=meth, cfg=cfg)
            @test sol_inplace.status == :Converged
            @test sol_inplace.residual_norm < res_atol
            @test maximum(abs.(sol_inplace.fixed_point .- legacy_simple_vec)) < val_atol
        end
    end

    # Parity on the Int -> Float promotion case
    parity_inputs_int = [1.0, 4.0]
    legacy_simple_vec2 = fixed_point(func, parity_inputs_int, LegacySimple()).fixed_point
    for (label, meth) in [("Simple_int_case", Simple()), ("Anderson_int_case", Anderson())]
        @testset "vector_parity_$label" begin
            sol = solve(func, copy(parity_inputs_int); method=meth, cfg=cfg)
            @test sol.status == :Converged
            @test sol.residual_norm < res_atol
            @test maximum(abs.(sol.fixed_point .- legacy_simple_vec2)) < val_atol
        end
    end
end

end
