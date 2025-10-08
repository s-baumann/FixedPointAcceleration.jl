using Random
using LinearAlgebra
using BenchmarkTools
using FixedPointAccelerationNext
using Printf

seed = 123
Random.seed!(seed)
N = 200
samples = 10

println("Legacy vs core comparison | N=$(N) samples=$(samples)")

function make_contraction(n::Int)
    d = 0.60 .+ 0.32 .* rand(n)
    A = Diagonal(d)
    b = randn(n)
    f(x) = A * x .+ b
    x0 = zeros(n)
    return f, x0
end

f, x0 = make_contraction(N)

cfg = FixedPointConfig(; threshold=1e-10, max_iters=5_000)

struct SolverStats
    iterations::Int
    residual::Float64
end

stats = Dict{Tuple{Symbol, Symbol}, SolverStats}()
median_times = Dict{Tuple{Symbol, Symbol}, Float64}()

using FixedPointAccelerationNext.OldImplementation:
    fixed_point,
    Simple as LegacySimple,
    Anderson as LegacyAnderson,
    Aitken as LegacyAitken,
    MPE as LegacyMPE,
    RRE as LegacyRRE

comparisons = [
    (:Simple, LegacySimple(), Simple()),
    (:Anderson, LegacyAnderson(; maxM=6), Anderson(; m=6)),
    (:Aitken, LegacyAitken(), Aitken()),
    (:MPE, LegacyMPE(; extrapolation_period=3), MPE(; period=3)),
    (:RRE, LegacyRRE(; extrapolation_period=3), RRE(; period=3)),
]

suite = BenchmarkGroup()

for (name, legacy_alg, core_method) in comparisons
    legacy_res = fixed_point(f, copy(x0), legacy_alg)
    stats[(:legacy, name)] = SolverStats(
        legacy_res.iterations,
        legacy_res.convergence_vector[end],
    )
    core_sol = solve(f, copy(x0); method=core_method, cfg=cfg)
    stats[(:core, name)] = SolverStats(core_sol.iterations, core_sol.residual_norm)

    group = BenchmarkGroup()
    group["legacy"] = @benchmarkable fixed_point($f, x, $legacy_alg) setup=(x = copy($x0);)
    group["core"] = @benchmarkable solve($f, x; method=$core_method, cfg=$cfg) setup=(x = copy($x0);)
    suite[string(name)] = group
end

results = BenchmarkTools.run(suite; samples=samples, evals=1, verbose=false)

for (name, _, _) in comparisons
    leg_trial = results[string(name)]["legacy"]
    core_trial = results[string(name)]["core"]
    median_times[(:legacy, name)] = BenchmarkTools.median(leg_trial).time / 1.0e9
    median_times[(:core, name)] = BenchmarkTools.median(core_trial).time / 1.0e9
end

println(
    rpad("Family", 8),
    rpad("Method", 10),
    rpad("Iter", 8),
    rpad("Time(s)", 12),
    "Residual",
)

for (name, _, _) in comparisons
    for family in (:legacy, :core)
        stat = stats[(family, name)]
        t = median_times[(family, name)]
        println(
            rpad(string(family), 8),
            rpad(string(name), 10),
            rpad(string(stat.iterations), 8),
            rpad(Printf.@sprintf("%.5f", t), 12),
            Printf.@sprintf("%.3e", stat.residual),
        )
    end
end

println("\nRelative iteration ratios (core / legacy):")
for (name, _, _) in comparisons
    leg = stats[(:legacy, name)]
    cor = stats[(:core, name)]
    ratio = cor.iterations / max(1, leg.iterations)
    println(rpad(string(name), 10), Printf.@sprintf("%.3f", ratio))
end

println("\nRelative median time ratios (core / legacy):")
for (name, _, _) in comparisons
    leg = median_times[(:legacy, name)]
    cor = median_times[(:core, name)]
    ratio = cor / max(1.0e-12, leg)
    println(rpad(string(name), 10), Printf.@sprintf("%.3f", ratio))
end

display(results)
