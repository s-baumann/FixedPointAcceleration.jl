module FixedPointAcceleration
using LinearAlgebra: cond, pinv
using Dates: now

include("options.jl")
include("results.jl")
include("function_execution.jl")
include("utilities.jl")

# Include all algorithm-specific files
include("algorithms/base.jl")
include("algorithms/simple.jl")
include("algorithms/anderson.jl")
include("algorithms/aitken.jl")
include("algorithms/newton.jl")
include("algorithms/mpe.jl")
include("algorithms/rre.jl")
include("algorithms/vea.jl")
include("algorithms/sea.jl")
include("algorithm_implementations.jl")  # This includes all the algorithm files

include("fixed_point.jl")

export fixed_point # Main functionality.

# Export algorithm types
export Simple, Anderson, Aitken, Newton, MPE, RRE, VEA, SEA

# Export configuration types and presets
export FixedPointOptions
export default_options, robust_options, fast_options, debug_options

end
