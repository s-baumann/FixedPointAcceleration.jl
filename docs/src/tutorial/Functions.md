# schumaker Function

The only exported function in this package (and the only one a basic user will ever need) is the schumaker function which creates a spline along with functions giving its derivative and second derivative. This is described as below

## Function documentation

```@docs
SchumakerSpline.schumaker
```

# Internal functions

## Imputing gradients

This section includes documentation on internal functions that are not exported from the package. A normal user will not need to look here but an advanced user may want to use internal functions directly or optimise the code for their particular project.

The first function imputes gradients if these are not provided.
```@docs
SchumakerSpline.imputeGradients
```

## Creating polynomials

These are used to create matrix that has the polynomial coefficients of the spline in each interval. Note that the extrapolate function is only used when linear or constant extrapolation is requested when outside the interpolation domain.

```@docs
SchumakerSpline.getCoefficientMatrix
SchumakerSpline.schumakerIndInterval
SchumakerSpline.extrapolate
```

## Creating splines from polynomials

Finally these functions create the spline function by executing the correct polynomial in each specific interval.

```@docs
SchumakerSpline.ppmak
SchumakerSpline.ppmakDeriv
SchumakerSpline.ppmak2Deriv
```
