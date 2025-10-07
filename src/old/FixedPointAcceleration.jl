# === Legacy implementation wrapped in submodule OldImplementation ===
module OldImplementation
using LinearAlgebra: cond, pinv
using GLM: LinearModel, fit
using Dates: now
include("options.jl")          # bring legacy option types into submodule scope
include("results.jl")
include("function_execution.jl")
include("utilities.jl")
# Legacy types (FixedPointOptions, FixedPointResults, presets) are defined locally
# via the includes above; no parent imports to avoid warnings.

"""Legacy implementation of FixedPointAcceleration preserved for backwards compatibility."""
# Include legacy algorithm stack isolated from top-level type names
include("algorithms/base.jl")
include("algorithms/simple.jl")
include("algorithms/anderson.jl")
include("algorithms/aitken.jl")
include("algorithms/newton.jl")
include("algorithms/mpe.jl")
include("algorithms/rre.jl")
include("algorithms/vea.jl")
include("algorithms/sea.jl")
include("algorithm_implementations.jl")
include("fixed_point.jl")

end # module OldImplementation

export OldImplementation
