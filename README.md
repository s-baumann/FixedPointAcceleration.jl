# FixedPointAcceleration.jl

| Build | Coverage | Documentation |
|-------|----------|---------------|
| [![Build Status](https://travis-ci.com/s-baumann/FixedPointAcceleration.jl.svg?branch=master)](https://travis-ci.org/s-baumann/FixedPointAcceleration.jl) | [![Coverage Status](https://coveralls.io/repos/github/s-baumann/FixedPointAcceleration.jl/badge.svg?branch=master)](https://coveralls.io/github/s-baumann/FixedPointAcceleration.jl?branch=master) | [![docs-latest-img](https://img.shields.io/badge/docs-latest-blue.svg)](https://s-baumann.github.io/FixedPointAcceleration.jl/dev/index.html) |

This package implements similar functionality to https://cran.r-project.org/web/packages/FixedPoint/index.html. The key differences are:
* This package makes use of Julia's type system and is generally typed in a more stable and extendible way. FixedPoint results are always output in a FixedPointResults struct. All algorithms are specified by an enum.
* This package does not include the plotting capability of the R package. This is not essential for the functionality of the package and since fixed_point calls can be chained together you can easily do whatever plotting you want pretty easily where necessary.

## Included Acceleration algorithms.

There are 8 acceleration algorithms included in this package. A more extensive explanation is available in the documentation for the R package:
https://cran.r-project.org/web/packages/FixedPoint/vignettes/FixedPoint.pdf
and the original papers are all cited there. A very brief description however is:
* Simple - This takes the output of the previous iterate and uses it as the next guess.

In addition the following three scalar algorithms can be used (they are elementwise for vectors):
* Aitken - This considers the sequence p, f(p), f(f(p)), ... convergences at a constant rate to the fixed point. After each two iterates it estimates the rate and jumps to the anticipated fixed point;
* Newton - This uses a Newtonian rootfinder to solve for the root of f(x) - x;
* SEA - or Scalar Epsilon Algorithm. This uses an epsilon algorithm for finding a fixed point where an elementwise inverse is taken;

In addition the following four algorithms are specialised for vectors of arbitrary length:
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

For more examples and information see the documentation.
