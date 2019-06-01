@enum FP_FunctionEvaluationError begin
    NoError = 0
    ErrorExecutingFunction = 1
    LengthOfOutputNotSameAsInput = 2
    InputMissingsDetected = 3
    InputNAsDetected = 4
    InputInfsDetected = 5
    OutputMissingsDetected = 6
    OutputNAsDetected = 7
    OutputInfsDetected = 8 # This is not always an error. The function f(x) = x has a fuxed point at infinity. In practical
    # settings however an infinite fixed point is unlikely to be the desired one. Also the algorithms here would not
    # likely work in the infinity case.
    FunctionIsNotTypeStable = 9
end

@enum FixedPointAccelerationAlgorithm begin
           Simple = 1
           Anderson = 2
           Newton = 3
           Aitken = 4
           MPE = 5
           RRE = 6
           SEA = 7
           VEA = 8
end

@enum InvalidReplacement begin
    NoAction = 1
    ReplaceElements = 2
    ReplaceVector = 3
end

@enum TerminationCondition begin
    ReachedConvergenceThreshold = 1
    ReachedMaxIter = 2
    InvalidInputOrOutputOfIteration = 3
end

struct FunctionEvaluationResult{RR<:Real,TT}
    Input_::Array{RR,1}
    Output_::Union{Missing,Array{Union{Missing,TT},1}}
    Other_Output_::Union{Missing,NamedTuple}
    Error_::FP_FunctionEvaluationError
    function FunctionEvaluationResult(Input::Array{T,1}, Output_::Missing, Error_::FP_FunctionEvaluationError, Other_Output_::Union{Missing,NamedTuple} = missing) where T<:Real
        return new{T,T}(Input, Output_, Other_Output_, Error_)
    end
    function FunctionEvaluationResult(Input::Array{T,1}, Output::Array{R,1}, Error_::FP_FunctionEvaluationError, Other_Output_::Union{Missing,NamedTuple} = missing) where T<:Real where R<:Real
        return new{T,R}(Input, convert(Array{Union{Missing,R},1}, Output), Other_Output_, Error_)
    end
    function FunctionEvaluationResult(Input::Array{T,1}, Output::Array{Union{Missing,R},1}, Error_::FP_FunctionEvaluationError, Other_Output_::Union{Missing,NamedTuple} = missing) where T<:Real where R<:Real
        return new{T,R}(Input, Output, Other_Output_, Error_)
    end
end

struct FixedPointResults{R<:Real}
    FixedPoint_::Union{Missing,Array{R,1}}
    Other_Output_::Union{Missing,NamedTuple}
    Convergence_::Union{Missing,R}          #
    TerminationCondition_::TerminationCondition
    Iterations_::Integer
    ConvergenceVector_::Union{Missing,Array{R,1}}
    FailedEvaluation_::Union{Missing,FunctionEvaluationResult}
    Inputs_::Array{R,2}
    Outputs_::Array{R,2}
    function FixedPointResults(Inputs_::Array{R,2}, Outputs_::Array{R,2}, TerminationCondition_::TerminationCondition;
                               ConvergenceVector_::Union{Missing,Array{R,1}} = missing,
                               FailedEvaluation_::Union{Missing,FunctionEvaluationResult} = missing,
                               Other_Output::Union{Missing,NamedTuple} = missing) where R<:Real
        Iterations_ = size(Outputs_)[2]
        FixedPoint_ = missing
        Convergence_ = missing
        if (!(ismissing(ConvergenceVector_))) && !(isempty(ConvergenceVector_))
            Convergence_ = ConvergenceVector_[Iterations_]
        end
        if TerminationCondition_ == ReachedConvergenceThreshold
            FixedPoint_ = Outputs_[:,Iterations_]
        end
        return new{R}(FixedPoint_, Other_Output, Convergence_, TerminationCondition_, Iterations_, ConvergenceVector_, FailedEvaluation_, Inputs_, Outputs_)
    end
end
