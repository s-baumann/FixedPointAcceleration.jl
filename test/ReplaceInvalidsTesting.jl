# Testing algos generating an invalid input.
function funcfunc(x::Array{Float64,1})
    # The first coordinate convergence to 4.0 by 1 unit per iterate.
    output = Array{Float64,1}(x)
    if abs(x[1] - 4.0) <= 1.0
        output[1] = 4.0
    elseif x[1] > 4.0
        output[1] = x[1] - 1.0
    else
        output[1] = x[1] + 1.0
    end
    # The second does aitken convergence to 2.3
    output[2] = x[2] + (2.3-x[2])/2.0
    return output
end
Inputs = [19.0,10.0]
fp = fixed_point(funcfunc, Inputs; Algorithm = Aitken)
fp.FailedEvaluation_.Error_ == InputInfsDetected
fp.FailedEvaluation_.Input_[1] == -Inf
!isinf(fp.FailedEvaluation_.Input_[2])
# Now fixing with replace element
fp = fixed_point(funcfunc, Inputs; Algorithm = Aitken, ReplaceInvalids = ReplaceElements)
fp.TerminationCondition_ == ReachedConvergenceThreshold
# Now fixing with replace element
fp2 = fixed_point(funcfunc, Inputs; Algorithm = Aitken, ReplaceInvalids = ReplaceVector)
fp2.TerminationCondition_ == ReachedConvergenceThreshold
fp.Iterations_ != fp2.Iterations_
