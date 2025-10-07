"""
 This contains the results from a function evaluation in the FixedPointAcceleration framework.
 It containing the following fields:
 * `Input_` - The input
 * `Output_` - The output of the  function. May be `missing` if function could not complete without error.
 * `Error_` - A symbol representing what error occured.
"""
struct FunctionEvaluationResult{T<:Number,R}
    input::Vector{T}
    output::Union{Missing,Vector{Missing},Vector{Union{Missing,R}},Vector{R}}
    other_output::Union{Missing,NamedTuple}
    error::Symbol
    function FunctionEvaluationResult(
        input_vec::Vector{T},
        output_val::Missing,
        error_sym::Symbol,
        other_output_val::Union{Missing,NamedTuple}=missing,
    ) where {T<:Number}
        return new{T,T}(input_vec, output_val, other_output_val, error_sym)
    end
    function FunctionEvaluationResult(
        input_vec::Vector{T},
        output_vec::Vector{R},
        error_sym::Symbol,
        other_output_val::Union{Missing,NamedTuple}=missing,
    ) where {T<:Number} where {R<:Union{Missing,<:Number}}
        if R === Missing
            return new{T,T}(input_vec, output_vec, other_output_val, error_sym)
        elseif R <: Number
            return new{T,R}(
                input_vec,
                convert(Vector{Union{Missing,R}}, output_vec),
                other_output_val,
                error_sym,
            )
        else
            return new{T,nonmissingtype(R)}(
                input_vec, output_vec, other_output_val, error_sym
            )
        end
    end
end

"""
 This contains the results from a the fixed_point function.

 It containing the following fields:
 * `FixedPoint_` - The `Vector` with the fixed point that has been found.
 * `Other_Output_` - The other output of the fixedpoint function.
 * `Convergence_` - A real number showing how close the `FixedPoint_` is to the input (to the function) that created it as output.
  * `TerminationCondition_' - Why did the fixedpoint acceleration stop.
  * `Iterations_' - How many iterations were undertaken
  * `ConvergenceVector_` - What is the convergence value at each iteration.
  * `FailedEvaluation_` - Why did the fixedpoint iteration fail (missing if it did not fail)
  * `Inputs_` - What were all of the inputs tried
  * `Outputs_` - What were all the corresponding outputs.
"""
struct FixedPointResults{R<:Number}
    fixed_point::Union{Missing,Vector{R}}
    other_output::Union{Missing,NamedTuple}
    convergence::Union{Missing,Real}
    termination_condition::Symbol
    iterations::Integer
    convergence_vector::Union{Missing,Vector{<:Real}}
    failed_evaluation::Union{Missing,FunctionEvaluationResult}
    inputs::Matrix{R}
    outputs::Matrix{R}
    function FixedPointResults(
        inputs_mat::Matrix{R},
        outputs_mat::Matrix{R},
        termination_condition::Symbol;
        convergence_vector::Union{Missing,Vector{<:Real}}=missing,
        failed_evaluation::Union{Missing,FunctionEvaluationResult}=missing,
        other_output_val::Union{Missing,NamedTuple}=missing,
    ) where {R<:Number}
        num_iterations = size(outputs_mat, 2)
        fixed_point_val = missing
        convergence_val = missing
        if !ismissing(convergence_vector) && !isempty(convergence_vector)
            convergence_val = convergence_vector[num_iterations]
        end
        if termination_condition == :ReachedConvergenceThreshold
            fixed_point_val = outputs_mat[:, num_iterations]
        end
        return new{R}(
            fixed_point_val,
            other_output_val,
            convergence_val,
            termination_condition,
            num_iterations,
            convergence_vector,
            failed_evaluation,
            inputs_mat,
            outputs_mat,
        )
    end
end
