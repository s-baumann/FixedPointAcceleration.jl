
# 3.  Using the FixedPoint package

## 3.1 Basic examples of using FixedPoint

### The Babylonian method for finding square roots.

Now we will demonstrate how **FixedPoint** can be used for simple problems. For the simplest possible case consider we want to estimate a square root using the Babylonian method. To find the square root of a number $x$, given an initial guess $t_0$ the following sequence converges to the square root:

$t_{n+1} = \frac{1}{2} \left[ t_n + \frac{x}{t_n} \right]$

This is a fast converging and inexpensive sequence which probably makes an acceleration algorithm overkill but for sake of exposition we can implement this in **FixedPoint**. In the next code block we find the square root of 100 with the simple method:

```
using FixedPointAcceleration
SequenceFunction = function(tn){0.5*(tn + 100/tn)}
Initial_Guess = 6.0
FP_Simple   = fixed_point(SequenceFunction,Initial_Guess; Method = "Simple")
```

We can also solve for a vector of fixed points at the same time. For instance every square root from 1 to 100.

```r
NumbersVector = 1:100
SequenceFunction = function(tn){0.5*(tn + NumbersVector/tn)}
FP_SEA   = FixedPoint(Function = SequenceFunction, Inputs = 1:100,  Method = "RRE")
```
Note that in this case the RRE method is being applied elementwise.

### Vectorised functions

The utility of the acceleration algorithms contained in **FixedPoint** are more apparent when applied to vectorised functions with cross dependency. For a simple example consider the below function where each entry of the vector depends on both entries of the previous iterate.

```r
SimpleVectorFunction = function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])}
FP_Simple = FixedPoint(Function = SimpleVectorFunction, Inputs = c(0.3,900), Method = "Simple")
FP_Anderson = FixedPoint(Function = SimpleVectorFunction, Inputs = c(0.3,900), Method = "Anderson")
```
This function takes 105 iterates to find a fixed point with the simple method but only 14 with the Anderson acceleration method.

## 3.2 Easily changing algorithm

A key priority in writing **FixedPoint** is for ease in changing fixed point method. The algorithm can be stopped and then restarted without having to repeat earlier iterations. This can be seen in the simple example below:

```r
Run1 = FixedPoint(Function = function(x){x^0.2}, Inputs = c(0.3, 0.4,0.5),
                Method = "Aitken")
cat(paste("Running this all at once takes ", length(Run1$Convergence), " iterations and
          convergences to ",  Run1$Convergence[length(Run1$Convergence)], "\n"))
```

```
## Running this all at once takes  9  iterations and
##           convergences to  1.11022302462516e-15
```

```r
Run2_1 = FixedPoint(Function = function(x){x^0.2}, Inputs = c(0.3, 0.4,0.5),
                Method = "Aitken", MaxIter = 4)
# Check some condition about the convergence or change algorithm.
Run2_2 = FixedPoint(Function = function(x){x^0.2}, Inputs = Run2_1$Inputs,
                    Outputs = Run2_1$Outputs, Method = "Aitken")
cat(paste("Running this with a break takes ", length(Run2_2$Convergence), " iterations and
          convergences to ",  Run2_2$Convergence[length(Run2_2$Convergence)], "\n"))
```

```
## Running this with a break takes  9  iterations and
##           convergences to  1.11022302462516e-15
```

This can be useful in problems where it is desired to switch extrapolation methods half way through. For instance in the consumption smoothing case presented later, the Vector Epsilon Algorithm converges quickly to a convergence level around $10^{-4}$ but then spends around 100 iterations hovering around this level without reducing further. In this case the algorithm could be changed to the simple method or Anderson method once the algorithm hit some target level of convergence.

Another example is parallelisation. The Anderson method in particular is well suited to use the results of previous iterates and these iterates do not need to be generated sequentially. If a suitable collection of guess vectors can be generated[^4] then all could be evaluated in parallel and recombined before again being used by the Anderson algorithm to generate a new guess. An example of this is presented later when the package is used to find the value function for a consumer in a consumption smoothing problem.

[^4]: For instance a set of guesses could be done by using different dampening values or by running the Anderson algorithm several times using different subsets of previous iterates to generate new guesses.

## 3.4 Graceful error handling

Hopefully **FixedPoint** is well tested enough that most kind of errors will be rare. **FixedPoint** also offers an option (``ReplaceInvalids'') to ensure that no acceleration algorithm generates guess vectors that contain NAs or Infs. This option can be set to ``ReplaceVector''  which will replace an extrapolated vector containing NAs or Infs by the vector output in the previous iterate. If it is set to ``ReplaceElement'' then it will replace the individual elements that are NA or Inf by the corresponding elements in the output of the previous iterate.[^15]

[^15]: With these options being set the **FixedPoint** should never throw an error when applied to a function which is guaranteed to map: $\Re^N \rightarrow \Re^N$ for some $N$ unless the error comes from underflow or overflow. If you do encounter such a testcase it likely reflects a bug, we would much appreciate it if you could send a minimal working example to the package maintainer.

Errors are likely however in cases where inputs functions have a restricted domain. For example this may include functions that require the input vector to have a particular shape (ie concavity) or functions where the input vector must be strictly positive. For a simple example consider the vectorised function we introduced in section 3.1. Now rather than $x^\prime[1] = \frac{\sqrt{\vert x[1] + x[2] \vert}}{2}$ we have $x^\prime[1] = \frac{\sqrt{ x[1] + x[2] }}{2}$ where the output x has a prime and the inputs has no prime symbol. $x^\prime[1]$ here is no longer real valued if the sum of the previous iterate is negative. This is what occurs in the 5th iterate of the Anderson method applied to this problem.

The FixedPoint function handles these situations gracefully by saving all previous results as well as the proposed new vector that lead to the error. In these cases the "Finish" message in the output describes the form of error, the "NewInputVector" list entry contains the attempted vector that was input to the function and the "NewOutputVector" list entry contains the result of applying the function to NewInputVector. This output information can be seen in the example below.

This information is useful in order to diagnose the issue. In this case we might decide to modify the function to insert the absolute value function with the reasoning that the same fixed point will likely result (which we could later verify). This also allows a user to run one method until an error occurs and then switch methods. This is demonstrated below.


```r
SimpleVectorFunction = function(x){c(0.5*sqrt(x[1] + x[2]), 1.5*x[1] + 0.5*x[2])}
FPSolution = FixedPoint(Function = SimpleVectorFunction, Inputs = c(0.3,900),
                        Method = "Anderson")
```

```
## Warning in sqrt(x[1] + x[2]): NaNs produced
```

```r
# We can see the output of the FixedPoint function in cases like this where it ends due
# to an error
FPSolution
```

```
## $Inputs
##       [,1]     [,2]      [,3]     [,4]
## [1,]   0.3  15.0025  7.350817 3.191655
## [2,] 900.0 450.4500 82.469248 9.574966
##
## $Outputs
##          [,1]      [,2]      [,3]     [,4]
## [1,]  15.0025  10.78717  4.738672 1.786520
## [2,] 450.4500 247.72875 52.260849 9.574966
##
## $Convergence
## [1] 449.550000 202.721250  30.208399   1.405135
##
## $FixedPoint
## [1] NA
##
## $Finish
## [1] "New output vector contains NAs"
##
## $NewInputVector
## [1] -1.085113 -3.255338
##
## $NewOutputVector
## [1]       NaN -3.255338
```

```r
# We can use this information to decide to switch to the simple method.
# No error results as the simple method doesn't extrapolate.
FPSolution = FixedPoint(Function = SimpleVectorFunction, Inputs = FPSolution$Inputs,
                        Outputs = FPSolution$Outputs, Method = "Simple", MaxIter = 5)
# Now we switch to the Anderson Method again. No error results because we are
# close to fixed point.
FPSolution = FixedPoint(Function = SimpleVectorFunction, Inputs = FPSolution$Inputs,
                        Outputs = FPSolution$Outputs, Method = "Anderson")
```



## 3.5 Convergence by constant increments

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
same set of fixed points. This is shown in the perceptron example in section 4.3. In other cases where it is not possible
to modify the function, it is advisable to use the simple method.

[^5]: When these complications arise the ReplaceInvalids method can be used to revert to a simple iterate or to change individual elements to the corresponding values in a simple iterate. This is as described in section 3.4.
