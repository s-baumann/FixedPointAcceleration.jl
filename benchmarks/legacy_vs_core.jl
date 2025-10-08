using Random
using LinearAlgebra
using BenchmarkTools
using FixedPointAccelerationNext
using FixedPointAcceleration
using Printf

seed = 123
Random.seed!(seed)
N = 200
samples = 10

println("Legacy vs new comparison | N=$(N) samples=$(samples)")

function make_contraction(n::Int)
    d = 0.60 .+ 0.32 .* rand(n)
    A = Diagonal(d) |> Matrix
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

stats = Dict{Tuple{Symbol,Symbol},SolverStats}()
median_times = Dict{Tuple{Symbol,Symbol},Float64}()

comparisons = [
    (:Simple, :Simple, Simple()),
    (:Anderson, :Anderson, Anderson(; m=6)),
    (:Aitken, :Aitken, Aitken()),
    (:MPE, :MPE, MPE(; period=3)),
    (:RRE, :RRE, RRE(; period=3)),
]

    legacy_res = fixed_point(
        f, copy(x0); Algorithm=:RRE, MaxM=6, ExtrapolationPeriod=3
    )


suite = BenchmarkGroup()

for (name, legacy_alg, new_method) in comparisons
    legacy_res = fixed_point(
        f, copy(x0); Algorithm=legacy_alg, MaxM=6, ExtrapolationPeriod=3
    )
    stats[(:legacy, name)] = SolverStats(
        legacy_res.Iterations_, legacy_res.ConvergenceVector_[end]
    )
    new_sol = solve(f, copy(x0); method=new_method, cfg=cfg)
    stats[(:new, name)] = SolverStats(new_sol.iterations, new_sol.residual_norm)

    group = BenchmarkGroup()
    group["legacy"] = @benchmarkable fixed_point($f, x; Algorithm=$legacy_alg, MaxM=6, ExtrapolationPeriod=3) setup = (x = copy($x0))
    group["new"] = @benchmarkable solve($f, x; method=$new_method, cfg=$cfg) setup = (
        x = copy($x0)
    )
    suite[string(name)] = group
end

results = BenchmarkTools.run(suite; samples, evals=5, verbose=true)

for (name, _, _) in comparisons
    leg_trial = results[string(name)]["legacy"]
    new_trial = results[string(name)]["new"]
    median_times[(:legacy, name)] = BenchmarkTools.median(leg_trial).time / 1.0e9
    median_times[(:new, name)] = BenchmarkTools.median(new_trial).time / 1.0e9
end

println(
    rpad("Family", 8), rpad("Method", 10), rpad("Iter", 8), rpad("Time(s)", 12), "Residual"
)

for (name, _, _) in comparisons
    for family in (:legacy, :new)
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

println("\nRelative iteration ratios (new / legacy):")
for (name, _, _) in comparisons
    leg = stats[(:legacy, name)]
    cor = stats[(:new, name)]
    ratio = cor.iterations / max(1, leg.iterations)
    println(rpad(string(name), 10), Printf.@sprintf("%.3f", ratio))
end

println("\nRelative median time ratios (new / legacy):")
for (name, _, _) in comparisons
    leg = median_times[(:legacy, name)]
    cor = median_times[(:new, name)]
    ratio = cor / max(1.0e-12, leg)
    println(rpad(string(name), 10), Printf.@sprintf("%.3f", ratio))
end

display(results)
