# FixedPointAcceleration.jl

| Build | Coverage |
|-------|----------|
| [![Build Status](https://travis-ci.com/s-baumann/FixedPointAcceleration.jl.svg?branch=master)](https://travis-ci.org/s-baumann/FixedPointAcceleration.jl) | [![Coverage Status](https://coveralls.io/repos/github/s-baumann/FixedPointAcceleration.jl/badge.svg?branch=master)](https://coveralls.io/github/s-baumann/FixedPointAcceleration.jl?branch=master)

This package implements similar functionality to https://cran.r-project.org/web/packages/FixedPoint/index.html. The key differences are:
* This package makes use of Julia's type system and is generally typed in a more stable and extendible way. FixedPoint results are always output in a FixedPointResults struct. All algorithms are specified by an enum.
* This package does not include the plotting capability of the R package. This is not essential for the functionality of the package and since fixed_point calls can be chained together you can easily do whatever plotting you want pretty easily where necessary.

## Included Acceleration algorithms.

There are 8 acceleration algorithms included in this package. A more extensive explanation is available in the documentation for the R package:
https://cran.r-project.org/web/packages/FixedPoint/vignettes/FixedPoint.pdf
and the original papers are all cited there. A very brief description however is:
* Simple - This takes the output of the previous iterate and uses it as the next guess.

In addition the following three scalar algorithms can be used (they are elementwise for vectors):
* Aitken - This considers the sequence p, f(p), f(f(p)) ... convergences at a constant rate to the fixed point. After each two iterates it estimates the rate and jumps to the anticipated fixed point;
* Newton - This uses a Newtonian rootfinder to solve for the root of f(x) - x;
* SEA - or Scalar Epsilon Algorithm. This uses an epsilon algorithm for finding a fixed point where an elementwise inverse is taken;

In addition the following four algorithms are specialised for vectors of arbitrary length. This is done after every "ExtrapolationPeriod" number of iterates which is 7 by default;
* VEA - or Vector Epsilon Algorithm. This uses an epsilon algorithm for finding a fixed point where an Moore-Penrose Pseudoinverse is used for vectors. This is done after every "ExtrapolationPeriod" number off iterates which is 7 by default;
* MPE - or Minimal Polynomial Extrapolation uses a linear combination of previous iterates for the next guess. This is done after every "ExtrapolationPeriod" number off iterates which is 7 by default;
* RRE - or Reduced Rank Extrapolation uses a linear combination of previous iterates for the next guess. This is done after every "ExtrapolationPeriod" number off iterates which is 7 by default;
* Anderson (default) - This takes a linear combination of previous iterates (it gets weights from an OLS regression). Unlike RRE, MPE,VEA this linear combination does not need to be of sequential simple iterates but can be any previous iterates. At maximum "MaxM" previous iterates are used but fewer may be used for numerical stability reasons if a certain matrix used in an OLS regression has a condition number above "ConditionNumberThreshold".


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
* All input functions must take and return Array{Float64,1} types. When solving for the fixedpoint of a scalar function use of the dot syntax - cos.(x) instead of cos(x) - might be useful for this.
* fp_anderson.FixedPoint_  in this case may be missing. This occurs whenever a fixedpoint was not found that is correct to the specified ConvergenceMetric and ConvergenceMetricThreshold. In this case the fp_anderson.Outputs_ array will probably contain something close to a fixed point in its rightmost column.

### Chaining Fixed Point Acceleration Algorithms

We can "chain" together different calls to the fixed_point function in order to switch acceleration algorithm
at any point. For instance consider the following function and initial guess at a fixed point:
```
func(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
Initial_Guess = [1.1,2.2]
```
Now we can initially do two simple iterates. Then do three iterates with the MPE method. Then one with the simple method and then finish with the RRE method. This can be done in the following way:
```
fp_chain      = fixed_point(func, Initial_Guess; Algorithm = Simple, MaxIter = 2)
fp_chain      = fixed_point(func, fp_chain; Algorithm = MPE, MaxIter = 3)
fp_chain      = fixed_point(func, fp_chain; Algorithm = Simple, MaxIter = 1)
fp_chain      = fixed_point(func, fp_chain; Algorithm = RRE, MaxIter = 100)
```
Now as it turns out The MPE (and RRE) does simple iterates except for every iterate that is a multiple of the ExtrapolationPeriod (7 by default). And so there is no difference from the above sequence of iterates and just doing all iterates with the RRE. This can be verified with the following:
```
fp_nochain = fixed_point(func, Inputs; Algorithm = RRE, MaxIter = 100)
fp_chain.Iterations_ == fp_nochain.Iterations_
all(abs.(fp_chain.Inputs_ .- fp_nochain.Inputs_) .< 1e-14)
```
This does highlight that there is no waste in changing fixed_point algorithm in this way. No iterates are reevaluated.

Changing algorithms can be useful in some cases where an error occurs. For instance consider we are trying to find the
fixed point of the following function:
```
simple_vector_function(x) = [0.5*sqrt(x[1] + x[2]), 1.5*x[1] + 0.5*x[2]]
Inputs = [0.3,900]
fp = fixed_point(simple_vector_function, Inputs; Algorithm = Anderson)
```
Inspecting this fp object reveals an error after the 3rditeration because Anderson tries to use a negative value for both x entries which results in the square root of a negative number. We can switch to simple iterations to get closer to the fixed point at which point Anderson will no longer try negative numbers. This will fix this.
```
fp = fixed_point(simple_vector_function, fp; Algorithm = Simple, MaxIter = 7)
fp = fixed_point(simple_vector_function, fp; Algorithm = Anderson)
```
## Termination conditions and Error handling.

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
