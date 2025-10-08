using BenchmarkTools
using FixedPointAccelerationNext

# Include individual benchmark files
include("simple_systems.jl")
include("optimal_growth_benchmark.jl")
include("poisson_benchmark.jl")
include("nash_equilibrium_benchmark.jl")
include("pagerank_benchmark.jl")
# include("legacy_vs_core.jl") # compare legacy vs new implementation

const SUITE = BenchmarkGroup()

simple_contraction!(SUITE)
nonlinear_system!(SUITE)
# optimal_growth!(SUITE)
# poisson!(SUITE)
# nash_equilibrium!(SUITE)
# pagerank!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
