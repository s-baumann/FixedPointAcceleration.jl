using Test
@testset "Test Chaining" begin
    using FixedPointAcceleration
    func1(x) = cos.(x)
    Inputs = 1.1
    # We should not have any slowdown from starting in simple for two. Then SEA for 3. Then Simple. Then VEA until the end.
    # The reason is that SEA/VEA does nothing of interest until 7 by default anyway.
    fp_chain = fixed_point(func1, Inputs, Simple(); MaxIter=Integer(2))
    fp_chain = fixed_point(func1, fp_chain, SEA(); MaxIter=Integer(3))
    fp_chain = fixed_point(func1, fp_chain, Simple(); MaxIter=Integer(1))
    fp_chain = fixed_point(func1, fp_chain, VEA(); MaxIter=Integer(100))

    fp_nochain = fixed_point(func1, Inputs, VEA(); MaxIter=Integer(100))
    @test fp_chain.Iterations_ == fp_nochain.Iterations_
    @test all(abs.(fp_nochain.Inputs_ .- fp_chain.Inputs_) .< 1e-14)

    func2(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
    Inputs = [1.1, 2.2]

    fp_chain = fixed_point(func2, Inputs, Simple(); MaxIter=Integer(2))
    fp_chain = fixed_point(func2, fp_chain, MPE(); MaxIter=Integer(3))
    fp_chain = fixed_point(func2, fp_chain, Simple(); MaxIter=Integer(1))
    fp_chain = fixed_point(func2, fp_chain, RRE(); MaxIter=Integer(100))

    fp_nochain = fixed_point(func2, Inputs, RRE(); MaxIter=Integer(100))
    @test fp_chain.Iterations_ == fp_nochain.Iterations_
    @test all(abs.(fp_nochain.Inputs_ .- fp_chain.Inputs_) .< 1e-14)
end
