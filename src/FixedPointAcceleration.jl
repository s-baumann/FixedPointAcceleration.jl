module FixedPointAcceleration
using LinearAlgebra: cond, pinv
using GLM: fit, LinearModel
using Dates: now

include("algorithms.jl")
include("types.jl")
include("function_execution.jl")
include("utilities.jl")
include("extrapolation_methods.jl")
include("algorithm_implementations.jl")
include("fixed_point_methods.jl")

export fixed_point, fixed_point_new_input # Main functionality.
export FunctionEvaluationResult, FixedPointResults # Structs.
export put_together_without_jumps, execute_function_safely # Auxillery functions. Mainly exported for testing.

# Export algorithm types
export FixedPointAlgorithm, Simple, Anderson, Aitken, Newton, MPE, RRE, VEA, SEA
export algorithm_name
export needs_extrapolation_period, get_extrapolation_period, is_polynomial_method, is_epsilon_method
end
