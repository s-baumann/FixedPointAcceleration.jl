struct FunctionEvaluationResult{RR<:Real,TT}
    Input_::Array{RR,1}
    Output_::Union{Missing,Array{Union{Missing,TT},1}}
    Other_Output_::Union{Missing,NamedTuple}
    Error_::Symbol
    function FunctionEvaluationResult(Input::Array{T,1}, Output_::Missing, Error_::Symbol, Other_Output_::Union{Missing,NamedTuple} = missing) where T<:Real
        return new{T,T}(Input, Output_, Other_Output_, Error_)
    end
    function FunctionEvaluationResult(Input::Array{T,1}, Output::Array{R,1}, Error_::Symbol, Other_Output_::Union{Missing,NamedTuple} = missing) where T<:Real where R<:Real
        return new{T,R}(Input, convert(Array{Union{Missing,R},1}, Output), Other_Output_, Error_)
    end
    function FunctionEvaluationResult(Input::Array{T,1}, Output::Array{Union{Missing,R},1}, Error_::Symbol, Other_Output_::Union{Missing,NamedTuple} = missing) where T<:Real where R<:Real
        return new{T,R}(Input, Output, Other_Output_, Error_)
    end
end

struct FixedPointResults{R<:Real}
    FixedPoint_::Union{Missing,Array{R,1}}
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
