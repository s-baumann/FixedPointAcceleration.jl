
# 3  Using the FixedPointAcceleration package

## 3.1 Basic examples of using FixedPointAcceleration

### The Babylonian method for finding square roots.

Now we will demonstrate how **FixedPointAcceleration** can be used for simple problems. For the simplest possible case consider we want to estimate a square root using the Babylonian method. To find the square root of a number $x$, given an initial guess $t_0$ the following sequence converges to the square root:

$t_{n+1} = \frac{1}{2} \left[ t_n + \frac{x}{t_n} \right]$

This is a fast converging and inexpensive sequence which probably makes an acceleration algorithm overkill but for sake of exposition we can implement this in **FixedPointAcceleration**. In the next code block we find the square root of 100 with the simple method:

```
using FixedPointAcceleration
SequenceFunction(x) = 0.5 .* (x .+ 100 ./ x)
Initial_Guess = 6.0
FP_Simple   = fixed_point(SequenceFunction, Initial_Guess; Algorithm = Simple)
```

We can also solve for a vector of fixed points at the same time. For instance every square root from 1 to 100.

```
NumbersVector = collect(1:100)
SequenceFunction(x) = 0.5 .* (x .+ NumbersVector ./ x)
Initial_Guess = repeat([10],100)
FP_SEA   = fixed_point(SequenceFunction, Initial_Guess; Algorithm = RRE)
```
Note that in this case the RRE method is being applied elementwise.

### Vectorised functions

The utility of the acceleration algorithms contained in **FixedPoint** are more apparent when applied to vectorised functions with cross dependency. For a simple example consider the below function where each entry of the vector depends on both entries of the previous iterate.

```
SimpleVectorFunction(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
Initial_Guess =  [0.3,900]
FP_Simple = fixed_point(SimpleVectorFunction  , Initial_Guess; Algorithm = Simple)
FP_Anderson = fixed_point(SimpleVectorFunction, Initial_Guess; Algorithm = Anderson)
```
This function takes 105 iterates to find a fixed point with the simple method but only 14 with the Anderson acceleration method.

## 3.2 Easily changing algorithm

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

## 3.3 Graceful error handling

Hopefully **FixedPointAcceleration** is well tested enough that most kind of errors will be rare. **FixedPointAcceleration** also offers an option (ReplaceInvalids) to ensure that no acceleration algorithm generates guess vectors that contain NaNs, Missings or Infs. This option can be set to ReplaceVector  which will replace an extrapolated vector containingNaNs, Missings or Infs by the vector output in the previous iterate. If it is set to ReplaceElement then it will replace the individual elements that are missings, NANs or Infs by the corresponding elements in the output of the previous iterate.

Errors are likely however in cases where inputs functions have a restricted domain. For example this may include functions that require the input vector to have a particular shape (ie concavity) or functions where the input vector must be strictly positive. For a simple example consider the vectorised function we introduced in section 3.1. Now rather than

$x^\prime[1] = \frac{\sqrt{\vert x[1] + x[2] \vert}}{2}$

we have

$x^\prime[1] = \frac{\sqrt{ x[1] + x[2] }}{2}$

where the output $$x$$ has a prime and the inputs has no prime symbol. $x^\prime[1]$ here is no longer real valued if the sum of the previous iterate is negative. This is what occurs in the 5th iterate of the Anderson method applied to this problem.

The FixedPoint function handles these situations gracefully by saving all previous results as well as the proposed new vector that lead to the error. In the event of such an error the FailedEvaluation_ member of the returned FixedPointResults struct will describe the issue.

This information is useful in order to diagnose the issue. In this case we might decide to modify the function to insert the absolute value function with the reasoning that the same fixed point will likely result (which we could later verify). This also allows a user to run one method until an error occurs and then switch methods. This is demonstrated below.

```
SimpleVectorFunction(x) = [0.5*sqrt(x[1] + x[2]), 1.5*x[1] + 0.5*x[2]]
Initial_Guess = [0.3,900]
FPSolution = FixedPoint(SimpleVectorFunction, Initial_Guess; Algorithm = Anderson)
```

```
# We can use this information to decide to switch to the simple method.
# No error results as the simple method doesn't extrapolate.
FPSolution = FixedPoint(SimpleVectorFunction, FPSolution; Algorithm = Simple, MaxIter = 5)
# Now we switch to the Anderson Method again. No error results because we are
# close to fixed point.
FPSolution = FixedPoint(SimpleVectorFunction, FPSolution; Algorithm = Anderson)
```

## 3.4 Convergence by constant increments

Most of the methods included in this function will fail in finding the fixed point of a function that converges by a fixed increment.
For instance we may have a function that takes $x$ and returns $x$ shifted 1 unit (in Euclidian norm) in a straight line
towards its fixed point. A realistic example of this is the training of a perceptron classifier which is explored later in section 4.3.

This case is problematic for all methods except for the simple method. The basic problem
can be illustrated simply by looking at the Newton method and the Aitken method. For the Newton method the derivative
is approximated by $\frac{ g(x_i) - g(x_{i-1})}{x_i-x{i-1}}$. When there is convergence by constant increments then
$g(x_i) = g(x_{i-1}) $ and the derivative is zero which means calculating the Newton method's recommended new guess of the
fixed point involves division by zero. Now considering the Aitken method the new guess is given by
$x_{i+1} = x_{i} - \frac{  (x_{i+1} - x_i)^2  }{  x_{i+2} - 2x_{i+1} + x_i}$.
When there is convergence by constant increments then $x_i - x_{i+1} = x_{i+1} - x_{i+2}$  and so we have $x_{i+2} - 2x_{i+1} + x_i = (x_i - x_{i+1}) - (x_{i+1} - x_{i+2}) = 0$. It is against not possible to calculate the new guess.[^5]

More generally similar problems exist for the other acceleration methods. When there is convergence by constant increments
then then the fixed point method receives information about what direction to go in but no information about how far to go.
This is a complication that is common to all of these acceleration methods in this package.
In these cases it may be possible to change the function to make it converge by varying increments while retaining the
same set of fixed points. This is shown in the perceptron example in section 4.2. In other cases where it is not possible
to modify the function, it is advisable to use the simple method.

[^5]: When these complications arise the ReplaceInvalids method can be used to revert to a simple iterate or to change individual elements to the corresponding values in a simple iterate. This is as described in section 3.3.
