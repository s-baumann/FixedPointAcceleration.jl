# Add standard contraction mapping benchmark
function simple_contraction!(suite)
    # Linear contraction mapping
    f(x) = 0.7 .* x .+ 1.2
    x0 = randn(50)
    cfg = FixedPointConfig(; threshold=1e-12, max_iters=1000)

    methods = [
        ("Simple", Simple()),
        ("Anderson", Anderson(; m=5)),
        ("Aitken", Aitken()),
        ("MPE", MPE(; period=5)),
        ("RRE", RRE(; period=5)),
        ("VEA", VEA(; period=7)),
        ("SEA", SEA(; period=7)),
    ]

    for (name, method) in methods
        suite["Polynomial system"]["Linear"][name] = @benchmarkable solve(
            $f, x0_copy; method=($method), cfg=($cfg)
        ) setup=(x0_copy = copy($x0))
    end
end

# Add nonlinear system benchmark
function nonlinear_system!(suite)
    # Nonlinear system: component-wise operations that create interdependencies
    function f!(dx, x::Vector{T}) where {T}
        n = length(x)
        for i in 1:n
            # Nonlinear coupling between components
            prev_idx = i == 1 ? n : i-1
            next_idx = i == n ? 1 : i+1
            dx[i] = 0.3 * x[i] + 0.1 * sin(x[prev_idx]) + 0.05 * x[next_idx]^2 + 0.2
        end
    end

    x0 = 0.5 * randn(30)
    cfg = FixedPointConfig(; threshold=1e-10, max_iters=1500)

    methods = [
        ("Simple", Simple()),
        ("Anderson", Anderson()),
        ("Aitken", Aitken()),
        ("MPE", MPE( )),
        ("RRE", RRE()),
        ("VEA", VEA()),
        ("SEA", SEA()),
    ]

    for (name, method) in methods
        suite["Polynomial system"]["Sinus-chain"][name] = @benchmarkable solve!(
            $f!, x0_copy; method=($method), cfg=($cfg)
        ) setup=(x0_copy = copy($x0))
    end
end
