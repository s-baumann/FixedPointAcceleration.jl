# FixedPointAcceleration.jl

This package implements similar functionality to https://cran.r-project.org/web/packages/FixedPoint/index.html. The key differences are:
* This package makes use of Julia's type system and is generally typed in a more stable way. FixedPoint results are always output in a FixedPointResults struct. All algorithms are specified by an enum.
* This package does not include the plots capability of the R package. This is not essential for the functionality of the package and since fixed_point calls can be chained together you can emulate regular plotting pretty easily where it is necessary.

## Using FixedPointAcceleration

We can first use FixedPointAcceleration to find the fixed point of the 1d function cos(x) with the Anderson method.
```
using FixedPointAcceleration
cos_func(x) = cos.(x)
Inputs = 1.1
fp_anderson = fixed_point(cos_func, Inputs; Algorithm = Anderson)
# And we can see the fixed point by looking at
fp_anderson.FixedPoint_
```
Two important issues are important to highlight here:
* All input functions must take and return Array{Float64,1} types. So when solving for the fixedpoint of a scalar function use the dot syntax - cos.(x) instead of cos(x). This is because the fixed_point function will try to evaluate the function with Array{Float64,1}([x]) rather than with x itself.
* fp_anderson.FixedPoint_  in this case may be missing. This occurs whenever a fixedpoint was not found that is correct to the specified ConvergenceMetric and ConvergenceMetricThreshold. In this case the fp_anderson.Outputs_ array will probably contain something close to a fixed point in its rightmost column.


Note that when

## Termination conditions and Error handling.

There are three possible TerminationCondition_ types that can be returned in the FixedPointResults struct. These are:
*  ReachedConvergenceThreshold - A fixed point has been reached.
*  ReachedMaxIter - The maximum number of iterations has been reached.
*  InvalidInputOrOutputOfIteration - A fatal error has occured while trying to solve for a fixed point. This is often simple to fix by simply changing algorithms for a while and hence any errors are caught and a FixedPointResults struct is returned detailing the error rather than explicitly throwing an error.

There are a few errors that can result in a InvalidInputOrOutputOfIteration termination. To aid in debugging where this termination condition is returned a FunctionEvaluationResult struct is returned as part of the FixedPointResults struct. This includes the inputs used when the error occured, the outputs (if they could be generated) and an additional error code (of enum type FP_FunctionEvaluationError):
* NoError - This indicates no error. You should never see this unless developing in the package as a function evaluation without an error will not cause a InvalidInputOrOutputOfIteration termination that causes the FunctionEvaluationResult struct to be returned.
* ErrorExecutingFunction - This indicates that there was an error evaluating the function with the given inputs. This will occur for instance if you try to evaluate sqrt.(x) at x = [-1.0] or 1/x at x = [0.0]. This may be solved by changing acceleration algorithm so that it does not try a vector which causes errors in the function. It may also be possible to reparameterise the function so that any vector is a valid input to the function.
* LengthOfOutputNotSameAsInput - A function taking an N-dimensional vector is not returning an N-dimensional vector.
* MissingsDetected - A function is returning an output vector containing missing values.
* NAsDetected - A function is returning an output vector containing NaN values.
* InfsDetected - A function is returning an output vector containing Inf values. While mathematically there is nothing wrong wtih this (Inf is a fixedpoint of the f(x) = x!), the algorithms of this package are not going to be useful in this case and hence it is not supported.

Together this error handling system should handle any errors gracefully without raising an ErrorException. ErrorExceptions are avoided so that the Inputs and Outputs from previous iterates are retained and the search for a fixed point can be resumed without interruption. If an ErrorException does occur while using fixed_point please raise an issue in github because this is not expected.
