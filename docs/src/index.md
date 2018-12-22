# SchumakerSpline

*A simple shape preserving spline implementation in Julia.*

A Julia package to create a shape preserving spline. This is a shape preserving spline which is guaranteed to be monotonic and concave/convex if the data is monotonic and concave/convex. It does not use any numerical optimisation and is therefore quick and smoothly converges to a fixed point in economic dynamics problems including value function iteration. It also automatically gives the first two derivatives
of the spline and options for determining behaviour when evaluated outside the interpolation domain.

This package has the same functionality as the R package called [schumaker](https://cran.r-project.org/web/packages/schumaker/index.html).

## Inputs

There are two optional setting in creating a spline. Firstly the gradients at each of the (x,y) points can be input to give more accuracy. If not supplied these are estimated from the points provided.

And secondly there are three options for out of sample prediction.

  * Curve - This is where the quadratic curve that is present in the first and last interval are used to predict points before the first interval and after the last interval respectively.

  * Linear - This is where a line is extended out before the first interval and after the last interval. The slope of the line is given by the derivative at the start of the first interval and end of the last interval.

  * Constant - This is where the first and last y values are used for prediction before the first point of the interval and after the last part of the interval respectively.


## Installation

```@contents
Pages = [
    "tutorials/Installation.md",
    ]
Depth = 2
```

## The functions of the package

```@contents
Pages = [
    "tutorials/Functions.md",
    ]
Depth = 2
```

## An example application

```@contents
Pages = [
    "tutorials/Example.md",
    ]
Depth = 2
```
