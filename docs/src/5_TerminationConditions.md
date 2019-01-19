# 5 Termination conditions and Error handling.

There are three possible TerminationCondition_ types that can be returned in the FixedPointResults struct. These are:
*  ReachedConvergenceThreshold - A fixed point has been reached.
*  ReachedMaxIter - The maximum number of iterations has been reached.
*  InvalidInputOrOutputOfIteration - A fatal error has occured while trying to solve for a fixed point. This is often simple to fix by simply changing algorithms for a while and hence any errors are caught and a FixedPointResults struct is returned detailing the error rather than explicitly throwing an error.

There are a few errors that can result in a InvalidInputOrOutputOfIteration termination. To aid in debugging where this termination condition is returned a FunctionEvaluationResult struct is returned as part of the FixedPointResults struct. This includes the inputs used when the error occured, the outputs (if they could be generated) and an additional error code (of enum type FP_FunctionEvaluationError):
* NoError - This indicates no error. You should never see this unless developing in the package as a function evaluation without an error will not cause a InvalidInputOrOutputOfIteration termination that causes the FunctionEvaluationResult struct to be returned.
* ErrorExecutingFunction - This indicates that there was an error evaluating the function with the given inputs. This will occur for instance if you try to evaluate sqrt.(x) at x = [-1.0] or 1/x at x = [0.0]. This may be solved by changing acceleration algorithm so that it does not try a vector which causes errors in the function. It may also be possible to reparameterise the function so that any vector is a valid input to the function.
* LengthOfOutputNotSameAsInput - A function taking an N-dimensional vector is not returning an N-dimensional vector.
* InputMissingsDetected - A function is returning an input vector containing missing values.
* InputNAsDetected - A function is returning an input vector containing NaN values.
* InputInfsDetected - A function is returning an input vector containing Inf values. While mathematically there is nothing wrong with this (Inf is a fixedpoint of the f(x) = x!), the algorithms of this package are not going to be useful in this case and hence it is not supported.
* OutputMissingsDetected - A function is returning an output vector containing missing values.
* OutputNAsDetected - A function is returning an output vector containing NaN values.
* OutputInfsDetected - A function is returning an output vector containing Inf values. While mathematically there is nothing wrong with this (like for InputInfsDetected) it is not supported.

Together this error handling system should handle any errors gracefully without raising an ErrorException. ErrorExceptions are avoided so that the Inputs and Outputs from previous iterates are retained and the search for a fixed point can be resumed without interruption. If an ErrorException does occur while using fixed_point please raise an issue in github because this is not expected.
