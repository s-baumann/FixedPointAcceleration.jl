using Test
using FixedPointAccelerationNext
using FixedPointAcceleration
# Core already loaded in runtests; remove duplicate using to avoid warnings

func(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
@testset "Simple Vector Function" begin
    x0 = [0.3, 900.0]
    fp_simple = solve(func, x0; method=Simple())
    @test fp_simple.residual_norm < 1e-10
    fp_anderson = solve(func, x0; method=Anderson())
    @test fp_anderson.residual_norm < 1e-10
    fp_aitken = solve(func, x0; method=Aitken())
    @test fp_aitken.residual_norm < 1e-10
    # fp_newton = solve(func, x0; method= Newton())
    # @test fp_newton.residual_norm < 1e-10
    fp_VEA = solve(func, x0; method=VEA())
    @test_broken fp_VEA.residual_norm < 1e-10
    fp_SEA = solve(func, x0; method=SEA())
    @test fp_SEA.residual_norm < 1e-10
    fp_MPE = solve(func, x0; method=MPE())
    @test fp_MPE.residual_norm < 1e-10
    fp_RRE = solve(func, x0; method=RRE())
    @test fp_RRE.residual_norm < 1e-10
end

@testset "new vs old parity" begin
    # Configuration for new core
    cfg = FixedPointConfig(;
        threshold=1e-10, max_iters=2_000, relaxation=1.0, relaxation_reference=:Simple
    )
    # Use original float input vector as baseline
    parity_inputs = [0.3, 900.0]
    legacy_simple_vec = fixed_point(func, parity_inputs; Algorithm=:Simple).FixedPoint_
    val_atol = 1e-8
    res_atol = 1e-10

    for (label, meth_legacy, meth_new) in [
        ("Simple_parity", :Simple, Simple()),
        ("Anderson_parity", :Anderson, Anderson()),
        ("Aitken_parity", :Aitken, Aitken()),
        ("MPE_parity", :MPE, MPE()),
        ("RRE_parity", :RRE, RRE()),
        ("VEA_parity", :VEA, VEA()),
        ("SEA_parity", :SEA, SEA()),
    ]
        @testset "vector_parity_$label" begin
            sol = solve(func, copy(parity_inputs); method=meth_new, cfg=cfg)
            @test sol.status == :Converged
            @test sol.residual_norm < res_atol
            @test maximum(abs.(sol.fixed_point .- legacy_simple_vec)) < val_atol
        end
        @testset "vector_parity_$(label)_both" begin
            sol = solve(func, copy(parity_inputs); method=meth_new, cfg=cfg)
            @test sol.status == :Converged
            @test sol.residual_norm < res_atol
            legacy = fixed_point(func, parity_inputs; Algorithm=meth_legacy).FixedPoint_
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

    # # Parity on the Int -> Float promotion case
    # parity_inputs_int = [1.0, 4.0]
    # legacy_simple_vec2 = fixed_point(func, parity_inputs_int; Algorithm=:Simple).FixedPoint_
    # for (label, meth) in [("Simple_int_case", Simple()), ("Anderson_int_case", Anderson())]
    #     @testset "vector_parity_$label" begin
    #         sol = solve(func, copy(parity_inputs_int); method=meth, cfg=cfg)
    #         @test sol.status == :Converged
    #         @test sol.residual_norm < res_atol
    #         @test maximum(abs.(sol.fixed_point .- legacy_simple_vec2)) < val_atol
    #     end
    # end
end
