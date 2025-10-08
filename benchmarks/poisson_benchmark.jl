using BenchmarkTools
using FixedPointAccelerationNext
using LinearAlgebra
using SparseArrays

"""
2D Poisson Equation: ∇²u = f(x,y) on [0,1]² with Dirichlet boundary conditions.
Discretized using finite differences, creating a large sparse linear system
that can be solved as a fixed point problem.
"""

struct PoissonProblem{T<:Real}
    nx::Int
    ny::Int
    h::T
    source::Function
    boundary::Function
end

function PoissonProblem(
    nx::Int,
    ny::Int;
    source=(x, y) -> -2π² * sin(π*x) * sin(π*y),  # Analytical solution: sin(πx)sin(πy)
    boundary=(x, y, side) -> 0.0,  # Homogeneous Dirichlet
)
    h = 1.0 / (nx + 1)
    return PoissonProblem{Float64}(nx, ny, h, source, boundary)
end

function poisson_operator(prob::PoissonProblem)
    nx, ny, h = prob.nx, prob.ny, prob.h

    # Create grid points (interior only)
    x_grid = [i * h for i in 1:nx]
    y_grid = [j * h for j in 1:ny]

    # Right-hand side vector
    b = zeros(nx * ny)
    for j in 1:ny
        for i in 1:nx
            idx = (j-1) * nx + i
            x, y = x_grid[i], y_grid[j]
            b[idx] = -h^2 * prob.source(x, y)
        end
    end

    # Boundary conditions
    for j in 1:ny
        for i in 1:nx
            idx = (j-1) * nx + i
            x, y = x_grid[i], y_grid[j]

            # Add boundary contributions
            if i == 1  # left boundary
                b[idx] += prob.boundary(0.0, y, :left)
            end
            if i == nx  # right boundary
                b[idx] += prob.boundary(1.0, y, :right)
            end
            if j == 1  # bottom boundary
                b[idx] += prob.boundary(x, 0.0, :bottom)
            end
            if j == ny  # top boundary
                b[idx] += prob.boundary(x, 1.0, :top)
            end
        end
    end

    # Jacobi iteration matrix
    return function (u::Vector{T}) where {T}
        u_new = similar(u)

        for j in 1:ny
            for i in 1:nx
                idx = (j-1) * nx + i

                # 5-point stencil
                sum_neighbors = 0.0

                # Left neighbor
                if i > 1
                    sum_neighbors += u[(j - 1) * nx + (i - 1)]
                end

                # Right neighbor
                if i < nx
                    sum_neighbors += u[(j - 1) * nx + (i + 1)]
                end

                # Bottom neighbor
                if j > 1
                    sum_neighbors += u[(j - 2) * nx + i]
                end

                # Top neighbor
                if j < ny
                    sum_neighbors += u[j * nx + i]
                end

                # Jacobi update
                u_new[idx] = 0.25 * (sum_neighbors + b[idx])
            end
        end

        return u_new
    end
end

function create_poisson_benchmark(grid_size::Int=20)
    prob = PoissonProblem(grid_size, grid_size)
    T = poisson_operator(prob)

    # Initial guess: random
    u0 = 0.1 * randn(grid_size * grid_size)

    cfg = FixedPointConfig(; threshold=1e-10, max_iters=2000)

    return T, u0, cfg, prob
end

function poisson!(suite)
    T, u0, cfg, prob = create_poisson_benchmark(15)  # 15x15 grid for benchmarking

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
        suite["Real-life problem"]["Poisson 2D PDE"][name] = @benchmarkable solve(
            $T, u0_copy; method=($method), cfg=($cfg)
        ) setup=(u0_copy = copy($u0))
    end
end
