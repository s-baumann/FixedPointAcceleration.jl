#!/usr/bin/env julia
# Comparison benchmark: legacy API vs unified core API across selected methods.
# Focuses on iteration counts and simple wall-clock (elapsed) times for a set of
# representative contraction problems.
# Usage:
#   julia --project=benchmarks benchmarks/legacy_vs_core2.jl
# Optional ENV:
#   FPA_COMPARE_DIM=200
#   FPA_COMPARE_REPEATS=5
#   FPA_COMPARE_SEED=123

using Random
using LinearAlgebra
using FixedPointAcceleration
using Printf

seed = parse(Int, get(ENV, "FPA_COMPARE_SEED", "123"));
Random.seed!(seed)
N = parse(Int, get(ENV, "FPA_COMPARE_DIM", "200"))
repeats = parse(Int, get(ENV, "FPA_COMPARE_REPEATS", "10"))

println("Legacy vs core comparison | N=$(N) repeats=$(repeats)")

# Build a diagonal linear contraction with mild difficulty (spectral radius ~0.92)
function make_contraction(n)
    d = 0.60 .+ 0.32 .* rand(n)  # in (0.60, 0.92)
    A = Diagonal(d)
    b = randn(n)
    f(x) = A * x .+ b
    x0 = zeros(n)
    return f, x0
end

f, x0 = make_contraction(N)

cfg = FixedPointConfig(; threshold=1e-10, max_iters=5_000)

struct Row
    family::Symbol      # :legacy or :core
    method::Symbol      # algorithm name
    iterations::Int
    time_s::Float64
    residual::Float64
end

rows = Row[]

# Helper to time a single solve (averaging repeats for smoother wall time)
function time_core(f, x0; method, cfg)
    best = Inf;
    sol_best = nothing
    for _ in 1:repeats
        t = @elapsed begin
            sol = solve(f, x0; method=method, cfg=cfg)
            sol_best = sol
        end
        best = t < best ? t : best
    end
    return sol_best, best
end

function time_legacy(f, x0; algorithm)
    best = Inf;
    sol_best = nothing
    for _ in 1:repeats
        t = @elapsed begin
            res = fixed_point(f, x0, algorithm)
            sol_best = res
        end
        best = t < best ? t : best
    end
    return sol_best, best
end

using FixedPointAcceleration.OldImplementation:
    fixed_point,
    Simple as LegacySimple,
    Anderson as LegacyAnderson,
    Aitken as LegacyAitken,
    MPE as LegacyMPE,
    RRE as LegacyRRE

# Legacy methods vs corresponding new core methods.
comparisons = [
    (:Simple, LegacySimple(), Simple()),
    (:Anderson, LegacyAnderson(; maxM=6), Anderson(; m=6)),
    (:Aitken, LegacyAitken(), Aitken()),
    (:MPE, LegacyMPE(; extrapolation_period=3), MPE(; period=3)),
    (:RRE, LegacyRRE(; extrapolation_period=3), RRE(; period=3)),
]

for (name, legacy_alg, core_method) in comparisons
    legacy_res, legacy_time = time_legacy(f, x0; algorithm=legacy_alg)
    push!(
        rows,
        Row(
            :legacy,
            name,
            legacy_res.iterations,
            legacy_time,
            legacy_res.convergence_vector[end],
        ),
    )
    core_sol, core_time = time_core(f, x0; method=core_method, cfg=cfg)
    push!(rows, Row(:core, name, core_sol.iterations, core_time, core_sol.residual_norm))
end

# Print summary
println(
    rpad("Family", 8), rpad("Method", 10), rpad("Iter", 8), rpad("Time(s)", 12), "Residual"
)
for r in rows
    println(
        rpad(string(r.family), 8),
        rpad(string(r.method), 10),
        rpad(string(r.iterations), 8),
        rpad(Printf.@sprintf("%.5f", r.time_s), 12),
        Printf.@sprintf("%.3e", r.residual)
    )
end

# Simple relative comparison
println("\nRelative iteration ratios (core / legacy):")
for name in (:Simple, :Anderson, :Aitken, :MPE, :RRE)
    leg = first(filter(row -> row.family == :legacy && row.method == name, rows))
    cor = first(filter(row -> row.family == :core && row.method == name, rows))
    ratio = cor.iterations / max(1, leg.iterations)
    println(rpad(string(name), 10), Printf.@sprintf("%.3f", ratio))
end
