module FixedPointAcceleration
    using LinearAlgebra: cond, pinv
    using GLM: fit, LinearModel
    using Dates: now

    include("types.jl")
    include("function_execution.jl")
    include("utilities.jl")
    include("extrapolation_methods.jl")
    include("algorithm_implementations.jl")
    include("fixed_point_methods.jl")
    export fixed_point, fixed_point_new_input # Main functionality.
    export FunctionEvaluationResult, FixedPointResults # Structs.
    export put_together_without_jumps, execute_function_safely # Auxillery functions. Mainly exported for testing.
end
