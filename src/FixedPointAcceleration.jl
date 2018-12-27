module FixedPointAcceleration
using LinearAlgebra
using GLM
include("0_Structs_And_Enums.jl")
include("MainFunctions.jl")
export fixed_point, fixed_point_new_input # Main functionality.
export FunctionEvaluationResult, FixedPointResults # Structs.
export put_together_without_jumps, create_safe_function_executor # Auxillery functions. Mainly exported for testing.
# Exporting enums
export FixedPointAccelerationAlgorithm, Simple, Anderson, Aitken, Newton, VEA, SEA, MPE, RRE # Algorithms
export InvalidReplacement, NoAction, ReplaceElements, ReplaceVector # Strategies for handling invalids in proposed new inputs.
export FP_FunctionEvaluationError, NoError, LengthOfOutputNotSameAsInput, MissingsDetected, NAsDetected, InfsDetected # Errors in evaluating
export TerminationCondition, AlreadyFixedPoint, ReachedConvergenceThreshold, ReachedMaxIter, InvalidInputOrOutputOfIteration # Termination Conditions
end
