"""
 This contains the results from a function evaluation in the FixedPointAcceleration framework.
 It containing the following fields:
 * `Input_` - The input
 * `Output_` - The output of the  function. May be `missing` if function could not complete without error.
 * `Error_` - A symbol representing what error occured.
"""
struct FunctionEvaluationResult{T<:Real,R}
    Input_::Vector{T}
    Output_::Union{Missing,Vector{Missing},Vector{Union{Missing,R}},Vector{R}}
    Other_Output_::Union{Missing,NamedTuple}
    Error_::Symbol
    function FunctionEvaluationResult(Input::Vector{T}, Output_::Missing, Error_::Symbol, Other_Output_::Union{Missing,NamedTuple} = missing) where T<:Real
        return new{T,T}(Input, Output_, Other_Output_, Error_)
    end
    function FunctionEvaluationResult(Input::Vector{T}, Output::Vector{R}, Error_::Symbol, Other_Output_::Union{Missing,NamedTuple} = missing) where T<:Real where R<:Union{Missing,<:Real}
        if R === Missing
            return new{T,T}(Input, Output, Other_Output_, Error_)
        elseif R <: Real
            return new{T,R}(Input, convert(Array{Union{Missing,R},1}, Output), Other_Output_, Error_)
        else
            return new{T,nonmissingtype(R)}(Input, Output, Other_Output_, Error_)
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
struct FixedPointResults{R<:Real}
    FixedPoint_::Union{Missing,Vector{R}}
    Other_Output_::Union{Missing,NamedTuple}
    Convergence_::Union{Missing,R}          #
    TerminationCondition_::Symbol
    Iterations_::Integer
    ConvergenceVector_::Union{Missing,Array{R,1}}
    FailedEvaluation_::Union{Missing,FunctionEvaluationResult}
    Inputs_::Array{R,2}
    Outputs_::Array{R,2}
    function FixedPointResults(Inputs_::Array{R,2}, Outputs_::Array{R,2}, TerminationCondition_::Symbol;
                               ConvergenceVector_::Union{Missing,Array{R,1}} = missing,
                               FailedEvaluation_::Union{Missing,FunctionEvaluationResult} = missing,
                               Other_Output::Union{Missing,NamedTuple} = missing) where R<:Real
        Iterations_ = size(Outputs_)[2]
        FixedPoint_ = missing
        Convergence_ = missing
        if (!(ismissing(ConvergenceVector_))) && !(isempty(ConvergenceVector_))
            Convergence_ = ConvergenceVector_[Iterations_]
        end
        if TerminationCondition_ == :ReachedConvergenceThreshold
            FixedPoint_ = Outputs_[:,Iterations_]
        end
        return new{R}(FixedPoint_, Other_Output, Convergence_, TerminationCondition_, Iterations_, ConvergenceVector_, FailedEvaluation_, Inputs_, Outputs_)
    end
end
