using FixedPointAcceleration
# Testing Error Evaluating Function
simple_vector_function(x) = [0.5*sqrt(x[1] + x[2]), 1.5*x[1] + 0.5*x[2]]
Inputs = [0.3,900]
fp = fixed_point(simple_vector_function, Inputs; Algorithm = Anderson)
# Inspecting this fp reveals an error after the 3rd iteration because
# Anderson tries to use a negative value for both x entries which results in
# the square root of a negative number. We can switch to simple
# iterations for a while to fix this.
fp.TerminationCondition_ == InvalidInputOrOutputOfIteration
fp.FailedEvaluation_.Error_ == ErrorExecutingFunction
fp = fixed_point(simple_vector_function, fp; Algorithm = Simple, MaxIter = Integer(7))
fp.TerminationCondition_ == ReachedMaxIter
fp = fixed_point(simple_vector_function, fp; Algorithm = Anderson)
fp.TerminationCondition_ == ReachedConvergenceThreshold

# Testing Input of NaN
fp = fixed_point(simple_vector_function, [NaN,900])
fp.FailedEvaluation_.Error_ == InputNAsDetected
# Testing Input of Inf
fp = fixed_point(simple_vector_function, [-Inf,900])
fp.FailedEvaluation_.Error_ == InputInfsDetected

# Testing Output of Nan
function funcfunc(x::Array{Float64,1})
    if abs(x[1] - 4.0) < 1e-12
        return Array{Float64,1}([NaN,4.0])
    end
    return sqrt.(x)
end
Inputs = [4.0,1.0]
fp = fixed_point(funcfunc, Inputs; Algorithm = Anderson)
fp.FailedEvaluation_.Error_ == OutputNAsDetected
# Testing Output of Missing
function funcfunc(x::Array{Float64,1})
    if abs(x[1] - 4.0) < 1e-12
        return [missing,4.0]
    end
    return sqrt.(x)
end
Inputs = [4.0,1.0]
fp = fixed_point(funcfunc, Inputs; Algorithm = Anderson)
fp.FailedEvaluation_.Error_ == OutputMissingsDetected
# Testing Output of Inf
function funcfunc(x::Array{Float64,1})
    if abs(x[1] - 4.0) < 1e-12
        return Array{Float64,1}([Inf,4.0])
    end
    return sqrt.(x)
end
Inputs = [4.0,1.0]
fp = fixed_point(funcfunc, Inputs; Algorithm = Anderson)
fp.FailedEvaluation_.Error_ == OutputInfsDetected
# Testing Output of wrong size
function funcfunc(x::Array{Float64,1})
    if abs(x[1] - 4.0) < 1e-12
        return Array{Float64,1}([5.0,4.0, 4.0])
    end
    return sqrt.(x)
end
Inputs = [4.0,1.0]
fp = fixed_point(funcfunc, Inputs; Algorithm = Anderson)
fp.FailedEvaluation_.Error_ == LengthOfOutputNotSameAsInput

# Testing Output of wrong type
function funcfunc(x::Array{Float64,1})
    return Array{Int,1}([5,4])
end
Inputs = [4.0,1.0]
fp = fixed_point(funcfunc, Inputs; Algorithm = Anderson)
fp.FailedEvaluation_.Error_ == FunctionIsNotTypeStable
