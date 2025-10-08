using BenchmarkTools
using FixedPointAccelerationNext

# include("legacy_vs_core.jl") # compare legacy vs new implementation

const SUITE = BenchmarkGroup()

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
