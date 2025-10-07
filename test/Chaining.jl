module test_Chaining

using Test
@testset "Test Chaining" begin
    using FixedPointAcceleration
    func1(x) = cos.(x)
    Inputs = 1.1
    # We should not have any slowdown from starting in simple for two. Then SEA for 3. Then Simple. Then VEA until the end.
    # The reason is that SEA/VEA does nothing of interest until 7 by default anyway.
    opts_2 = FixedPointOptions(max_iterations=2)
    opts_3 = FixedPointOptions(max_iterations=3)
    opts_1 = FixedPointOptions(max_iterations=1)
    opts_100 = FixedPointOptions(max_iterations=100)

    fp_chain = fixed_point(func1, Inputs, Simple(), opts_2)
    fp_chain = fixed_point(func1, fp_chain, SEA(), opts_3)
    fp_chain = fixed_point(func1, fp_chain, Simple(), opts_1)
    fp_chain = fixed_point(func1, fp_chain, VEA(), opts_100)

    fp_nochain = fixed_point(func1, Inputs, VEA(), opts_100)
    @test fp_chain.Iterations_ == fp_nochain.Iterations_
    @test all(abs.(fp_nochain.Inputs_ .- fp_chain.Inputs_) .< 1e-14)

    func2(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
    Inputs = [1.1, 2.2]

    fp_chain = fixed_point(func2, Inputs, Simple(), opts_2)
    fp_chain = fixed_point(func2, fp_chain, MPE(), opts_3)
    fp_chain = fixed_point(func2, fp_chain, Simple(), opts_1)
    fp_chain = fixed_point(func2, fp_chain, RRE(), opts_100)

    fp_nochain = fixed_point(func2, Inputs, RRE(), opts_100)
    @test fp_chain.Iterations_ == fp_nochain.Iterations_
    @test all(abs.(fp_nochain.Inputs_ .- fp_chain.Inputs_) .< 1e-14)
end
end
