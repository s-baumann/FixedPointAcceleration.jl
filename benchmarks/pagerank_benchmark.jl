using BenchmarkTools
using FixedPointAccelerationNext
using LinearAlgebra
using SparseArrays
using Random

"""
PageRank Algorithm: Find the stationary distribution of a Markov chain.
This creates a large-scale sparse linear algebra problem that's excellent
for testing acceleration methods on practical applications.
"""

function create_random_graph(n::Int; density::Float64=0.1, seed::Int=42)
    Random.seed!(seed)

    # Create random adjacency matrix
    A = spzeros(n, n)
    for i in 1:n
        for j in 1:n
            if i != j && rand() < density
                A[i, j] = rand()
            end
        end
    end

    # Ensure every node has at least one outgoing edge (add self-loops if needed)
    for i in 1:n
        if sum(A[i, :]) == 0
            A[i, i] = 1.0
        end
    end

    return A
end

function pagerank_operator(A::AbstractMatrix, α::Float64=0.85)
    n = size(A, 1)

    # Create transition matrix (column stochastic)
    P = similar(A)
    for j in 1:n
        col_sum = sum(A[:, j])
        if col_sum > 0
            P[:, j] = A[:, j] / col_sum
        else
            P[:, j] = fill(1/n, n)  # Uniform distribution for dangling nodes
        end
    end

    # PageRank operator: αP + (1-α)/n * ee^T
    teleport = (1 - α) / n

    return function (x::Vector{T}) where {T}
        # PageRank iteration: x_{k+1} = αPx_k + (1-α)/n * e
        return α * (P * x) .+ teleport
    end
end

function create_pagerank_benchmark(n::Int=100, density::Float64=0.08)
    A = create_random_graph(n; density=density)
    T = pagerank_operator(A, 0.85)  # Standard PageRank damping factor

    # Initial guess: uniform distribution
    x0 = fill(1.0/n, n)

    cfg = FixedPointConfig(; threshold=1e-12, max_iters=2000)

    return T, x0, cfg, A
end

function pagerank!(suite)
    T, x0, cfg, A = create_pagerank_benchmark(80, 0.1)  # 80 nodes, 10% edge density

    methods = [
        ("Simple", Simple()),
        ("Anderson", Anderson(; m=8)),
        ("Aitken", Aitken()),
        ("MPE", MPE(; period=6)),
        ("RRE", RRE(; period=6)),
        ("VEA", VEA(; period=8)),
        ("SEA", SEA(; period=8)),
    ]

    for (name, method) in methods
        suite["Real-life problem"]["Pagerank"][name] = @benchmarkable solve(
            $T, x0_copy; method=($method), cfg=($cfg)
        ) setup=(x0_copy = copy($x0))
    end
end
