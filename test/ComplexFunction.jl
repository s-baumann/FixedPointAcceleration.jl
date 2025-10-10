using Test

@testset "Complex Fixed Point Support" begin
    using FixedPointAcceleration

    func_scalar(z) = (z .+ (1 + 2im)) ./ 2
    scalar_start = 0.0 + 0.0im

    res_simple = fixed_point(func_scalar, scalar_start; Algorithm = :Simple, MaxIter = 200)
    @test res_simple.TerminationCondition_ == :ReachedConvergenceThreshold
    @test res_simple.Convergence_ < 1e-10
    @test res_simple.FixedPoint_ isa Vector{ComplexF64}
    @test isapprox(res_simple.FixedPoint_[1], 1 + 2im; atol = 1e-8)

    res_anderson = fixed_point(func_scalar, scalar_start; Algorithm = :Anderson, MaxIter = 200)
    @test res_anderson.TerminationCondition_ == :ReachedConvergenceThreshold
    @test isapprox(res_anderson.FixedPoint_[1], 1 + 2im; atol = 1e-8)

    func_vec(v) = (v .+ ComplexF64[1 + 2im, -3 + 1im]) ./ 2
    vec_start = ComplexF64[0.25 - 0.5im, 0.75 + 0.1im]
    res_vec = fixed_point(func_vec, vec_start; Algorithm = :Anderson, MaxIter = 200)
    @test res_vec.TerminationCondition_ == :ReachedConvergenceThreshold
    @test all(isapprox.(res_vec.FixedPoint_, ComplexF64[1 + 2im, -3 + 1im]; atol = 1e-8))

    new_guess = fixed_point_new_input(res_vec.Inputs_, res_vec.Outputs_, :Anderson)
    @test eltype(new_guess) <: Complex
    @test !any(isnan.(new_guess))

    resumed = fixed_point(func_vec, res_vec; Algorithm = :Anderson)
    @test resumed.TerminationCondition_ == :ReachedConvergenceThreshold
    @test all(isapprox.(resumed.FixedPoint_, ComplexF64[1 + 2im, -3 + 1im]; atol = 1e-8))
end
