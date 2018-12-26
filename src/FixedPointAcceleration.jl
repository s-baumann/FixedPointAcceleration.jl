module FixedPointAcceleration
using LinearAlgebra
using GLM
include("FixedPointResults.jl")
include("MainFunctions.jl")
export fixed_point, fixed_point_new_input
export put_together_without_jumps, create_safe_function_executor
end
