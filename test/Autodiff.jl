using Test
using ForwardDiff
using FixedPointAcceleration
using FixedPointAcceleration: solve, FixedPointConfig, Simple

function model_func(
    x::AbstractVector{T}, params::AbstractVector{T}
) where {T<:Real}
    return [params[1] * sqrt(abs(x[1] + x[2])), 1.5 * x[1] + params[2] * x[2]]
end

function ML_estimation(params::AbstractVector{T}) where {T<:Real}
    inputs = T[zero(T) + 0.3, zero(T) + 900.0]
    sol = solve(
        inp -> model_func(inp, params),
        inputs;
        method=Simple(),
        cfg=FixedPointConfig(; threshold=1e-12),
    )
    return sol.fixed_point[1]
end

@testset "Automatic Differentiation" begin
    params = [1.5, 0.5]
    grads = ForwardDiff.gradient(ML_estimation, params)
    @test isa(grads, Vector{Float64})
    @test all(isfinite, grads)
end
