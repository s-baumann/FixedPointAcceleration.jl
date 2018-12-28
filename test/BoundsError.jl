# Testing Error Evaluating Function
simple_vector_function(x) = [0.5*sqrt(x[1] + x[2]), 1.5*x[1] + 0.5*x[2]]
Inputs = [0.3,900]
fp = fixed_point(simple_vector_function, Inputs; Algorithm = Anderson)
# Inspecting this fp reveals an error after the 3rditeration because
# Anderson tries to use a negative value for both x entries which results in
# the square root of a negative number. We can switch to simple
# iterations for a while to fix this.
fp.TerminationCondition_ == InvalidInputOrOutputOfIteration
fp.FailedEvaluation_.Error_ == ErrorExecutingFunction
fp
fp = fixed_point(simple_vector_function, fp; Algorithm = Simple, MaxIter = 7)
fp.TerminationCondition_ == ReachedMaxIter
fp = fixed_point(simple_vector_function, fp; Algorithm = Anderson)
fp.TerminationCondition_ == ReachedConvergenceThreshold

# Testing Output of Nan
#function funcfunc(x::Array{Float64,1})
#    if abs(x[1] - 4.0) < 1e-12
#        return [NaN,NaN]
#    end
#    return sqrt.(x)
#end
#Inputs = [4.0,1.0]
#fp = fixed_point(funcfunc, Inputs; Algorithm = Anderson)
#fp.FailedEvaluation_.Error_
