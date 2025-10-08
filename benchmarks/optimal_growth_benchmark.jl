using BenchmarkTools
using FixedPointAccelerationNext
using LinearAlgebra
using Interpolations

"""
Optimal Growth Model with CRRA utility and Cobb-Douglas production.
This creates a high-dimensional fixed point problem that tests various acceleration methods.
"""

struct OptimalGrowthModel{T<:Real}
    β::T        # discount factor
    σ::T        # risk aversion
    α::T        # capital share
    δ::T        # depreciation rate
    k_grid::Vector{T}  # capital grid
    n_k::Int    # grid size
end

function OptimalGrowthModel(; β=0.96, σ=2.0, α=0.36, δ=0.1, k_min=0.1, k_max=10.0, n_k=100)
    k_grid = collect(range(k_min, k_max; length=n_k))
    return OptimalGrowthModel(β, σ, α, δ, k_grid, n_k)
end

function bellman_operator(model::OptimalGrowthModel)
    return function (V::Vector{T}) where {T}
        V_new = similar(V)

        for (i, k) in enumerate(model.k_grid)
            # Production and available resources
            y = k^model.α
            max_c = y + (1 - model.δ) * k

            # Find optimal consumption/saving
            max_value = -Inf

            for (j, k_next) in enumerate(model.k_grid)
                c = max_c - k_next
                if c > 0
                    # CRRA utility
                    u = c^(1 - model.σ) / (1 - model.σ)
                    # Expected continuation value (deterministic case)
                    continuation = model.β * V[j]
                    value = u + continuation

                    if value > max_value
                        max_value = value
                    end
                end
            end

            V_new[i] = max_value
        end

        return V_new
    end
end

function create_optimal_growth_benchmark()
    model = OptimalGrowthModel(; n_k=50)  # Smaller for benchmarking
    T = bellman_operator(model)

    # Initial guess: linear in capital
    V0 = log.(model.k_grid)

    cfg = FixedPointConfig(; threshold=1e-8, max_iters=1000)

    return T, V0, cfg, model
end

# Benchmark suite
function optimal_growth!(suite)
    T, V0, cfg, model = create_optimal_growth_benchmark()

    methods = [
        ("Simple", Simple()),
        ("Anderson", Anderson()),
        ("Aitken", Aitken()),
        ("MPE", MPE()),
        ("RRE", RRE()),
        ("VEA", VEA()),
        ("SEA", SEA()),
    ]

    for (name, method) in methods
        suite["Real-life problem"]["Optimal Growth"][name] = @benchmarkable solve(
            $T, V0_copy; method=($method), cfg=($cfg)
        ) setup=(V0_copy = copy($V0))
    end
end
