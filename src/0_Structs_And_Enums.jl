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

struct FunctionEvaluationResult
    Input_::Array{Float64,1}
    Output_::Union{Missing,Array{Float64,1}}
    Error_::FP_FunctionEvaluationError
    function FunctionEvaluationResult(Input::Array, Output_::Missing, Error_::FP_FunctionEvaluationError)
        return new(vec(Input), Output_, Error_)
    end
    function FunctionEvaluationResult(Input::Array, Output::Array, Error_::FP_FunctionEvaluationError)
        return new(vec(Input), vec(Output), Error_)
    end
end

struct FixedPointResults
    FixedPoint_::Union{Missing,Array{Float64,1}}
    Convergence_::Union{Missing,Float64}          #
    TerminationCondition_::TerminationCondition
    Iterations_::Int
    ConvergenceVector_::Union{Missing,Array{Float64,1}}
    FailedEvaluation_::Union{Missing,FunctionEvaluationResult}
    Inputs_::Array{Float64,2}
    Outputs_::Array{Float64,2}
    function FixedPointResults(Inputs_::Array{Float64,2}, Outputs_::Array{Float64,2}, TerminationCondition_::TerminationCondition;
                               ConvergenceVector_::Union{Missing,Array{Float64,1}} = missing,
                               FailedEvaluation_::Union{Missing,FunctionEvaluationResult} = missing)
        Iterations_ = size(Outputs_)[2]
        FixedPoint_ = missing
        Convergence_ = missing
        if (!(ismissing(ConvergenceVector_))) && !(isempty(ConvergenceVector_))
            Convergence_ = ConvergenceVector_[Iterations_]
        end
        if TerminationCondition_ == ReachedConvergenceThreshold
            FixedPoint_ = Outputs_[:,Iterations_]
        end
        return new(FixedPoint_, Convergence_, TerminationCondition_, Iterations_, ConvergenceVector_, FailedEvaluation_, Inputs_, Outputs_)
    end
end
