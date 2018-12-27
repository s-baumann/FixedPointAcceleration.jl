@enum FP_FunctionEvaluationError begin
    NoError = 0
    LengthOfOutputNotSameAsInput = 1
    MissingsDetected = 2
    NAsDetected = 3
    InfsDetected = 4 # This is not always an error. The function f(x) = x has a fuxed point at infinity. In practical
    # settings however an infinite fixed point is unlikely to be the desired one. Also the algorithms here would not
    # likely work in the infinity case.
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
    AlreadyFixedPoint = 0
    ReachedConvergenceThreshold = 1
    ReachedMaxIter = 2
    InvalidInputOrOutputOfIteration = 3
end

struct FunctionEvaluationResult
    AbortFunction::Bool
    Input::Array{Float64,1}
    Output::Array{Float64,1}
    Message::FP_FunctionEvaluationError
end

struct FixedPointResults
    FixedPoint_::Union{Missing,Array{Float64,1}}
    Convergence_::Union{Missing,Float64}          #
    TerminationCondition_::TerminationCondition
    Iterations_::Int
    ConvergenceVector_::Union{Missing,Array{Float64,1}}
    FailedEvaluation_::Union{Missing,NamedTuple}
    Inputs_::Array{Float64,2}
    Outputs_::Array{Float64,2}
    function FixedPointResults(Inputs_::Array{Float64,2}, Outputs_::Array{Float64,2}, TerminationCondition_::TerminationCondition;
                               ConvergenceVector_::Union{Missing,Array{Float64,1}} = missing,
                               FailedEvaluation_::Union{Missing,NamedTuple} = missing)
        Iterations_ = size(Inputs_)[2]
        FixedPoint_ = missing
        Convergence_ = missing
        if !(isempty(Outputs_))
            FixedPoint_ = Outputs_[:,Iterations_]
        end
        if (!(ismissing(ConvergenceVector_))) && !(isempty(ConvergenceVector_))
            Convergence_ = ConvergenceVector_[Iterations_]
        end
        return new(FixedPoint_, Convergence_, TerminationCondition_, Iterations_, ConvergenceVector_, FailedEvaluation_, Inputs_, Outputs_)
    end
end
