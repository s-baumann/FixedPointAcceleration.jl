module FixedPointAcceleration
using LinearAlgebra: cond, pinv
using GLM: fit, LinearModel
using Dates: now
include("0_Structs_And_Enums.jl")
include("1_MainFunctions.jl")
export fixed_point, fixed_point_new_input # Main functionality.
export FunctionEvaluationResult, FixedPointResults # Structs.
export put_together_without_jumps, execute_function_safely # Auxillery functions. Mainly exported for testing.
end
