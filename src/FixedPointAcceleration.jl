module FixedPointAcceleration
using LinearAlgebra
using GLM
include("FixedPointResults.jl")
include("MainFunctions.jl")
export fixed_point, fixed_point_new_input
export put_together_without_jumps, create_safe_function_executor
# Exporting enums
export FixedPointAccelerationAlgorithm, Simple, Anderson, Aitken, Newton, VEA, SEA, MPE, RRE
export InvalidReplacement, NoAction, ReplaceElements, ReplaceVector
export FixedPointAccelerationPlots, NoPlot, ConvergenceFig, ChangePerIterate

end
