
<a id='schumaker-Function-1'></a>

# schumaker Function


The only exported function in this package (and the only one a basic user will ever need) is the schumaker function which creates a spline along with functions giving its derivative and second derivative. This is described as below


<a id='Function-documentation-1'></a>

## Function documentation

<a id='SchumakerSpline.schumaker' href='#SchumakerSpline.schumaker'>#</a>
**`SchumakerSpline.schumaker`** &mdash; *Function*.



Creates splines for a given set of x and y values (and optionally gradients) and the first and second derivatives of this spline.

**Takes**

  * x - A Float64 vector of x coordinates.
  * y - A Float64 vector of y coordinates.
  * gradients (optional)- A Float64 vector of gradients at each point. If not supplied these are imputed from x and y.
  * extrapolation (optional) - This should be a string in ("Curve", "Linear", "Constant") specifying how to interpolate outside of the sample domain. By default it is "curve" which extends out the first and last quadratic curves. The other options are "Linear" which extends the line (from first and last curve) out from the first and last point and "Constant" which extends out the y value at the first and last point.

**Returns**

  * A spline which takes and input value and returns the spline y value.
  * The derivative of this spline.
  * The second derivative of this spline.  


<a id='Internal-functions-1'></a>

# Internal functions


<a id='Imputing-gradients-1'></a>

## Imputing gradients


This section includes documentation on internal functions that are not exported from the package. A normal user will not need to look here but an advanced user may want to use internal functions directly or optimise the code for their particular project.


The first function imputes gradients if these are not provided.

<a id='SchumakerSpline.imputeGradients' href='#SchumakerSpline.imputeGradients'>#</a>
**`SchumakerSpline.imputeGradients`** &mdash; *Function*.



Imputes gradients based on a vector of x and y coordinates.

**Takes**

  * x - A Float64 vector of x coordinates
  * y - A Float64 vector of y coordinates

**Returns**

  * A Float64 vector of gradients for each input point


<a id='Creating-polynomials-1'></a>

## Creating polynomials


These are used to create matrix that has the polynomial coefficients of the spline in each interval. Note that the extrapolate function is only used when linear or constant extrapolation is requested when outside the interpolation domain.

<a id='SchumakerSpline.getCoefficientMatrix' href='#SchumakerSpline.getCoefficientMatrix'>#</a>
**`SchumakerSpline.getCoefficientMatrix`** &mdash; *Function*.



Calls SchumakerIndInterval many times to get full set of spline intervals and coefficients. Then calls extrapolation for out of sample behaviour

**Takes**

  * gradients - A Float64 vector of gradients at each point
  * x - A Float64 vector of x coordinates
  * y - A Float64 vector of y coordinates
  * extrapolation - A string in ("Curve", "Linear", "Constant") that gives behaviour outside of interpolation range.

**Returns**

  * A vector of interval starts
  * A vector of interval ends
  * A matrix of all coefficients   

<a id='SchumakerSpline.schumakerIndInterval' href='#SchumakerSpline.schumakerIndInterval'>#</a>
**`SchumakerSpline.schumakerIndInterval`** &mdash; *Function*.



Splits an interval into 2 subintervals and creates the quadratic coefficients

**Takes**

  * s - A 2 entry Float64 vector with gradients at either end of the interval
  * z - A 2 entry Float64 vector with y values at either end of the interval
  * Smallt - A 2 entry Float64 vector with x values at either end of the interval

**Returns**

  * A 2 x 5 matrix. The first column is the x values of start of the two subintervals. The second column is the ends. The last 3 columns are quadratic coefficients in two subintervals.

<a id='SchumakerSpline.extrapolate' href='#SchumakerSpline.extrapolate'>#</a>
**`SchumakerSpline.extrapolate`** &mdash; *Function*.



Adds a row on top and bottom of coefficient matrix to give out of sample prediction.

**Takes**

  * fullMatrix - output from GetCoefficientMatrix first few lines
  * extrapolation - A string in ("Curve", "Linear", "Constant") that gives behaviour outside of interpolation range.
  * x - A Float64 vector of x coordinates
  * y - A Float64 vector of y coordinates

**Returns**

  * A new version of fullMatrix with out of sample prediction built into it.


<a id='Creating-splines-from-polynomials-1'></a>

## Creating splines from polynomials


Finally these functions create the spline function by executing the correct polynomial in each specific interval.

<a id='SchumakerSpline.ppmak' href='#SchumakerSpline.ppmak'>#</a>
**`SchumakerSpline.ppmak`** &mdash; *Function*.



Creates a spline defined by interval starts IntStarts and quadratic coefficients SpCoefs which evaluates an input point.

**Takes**

  * IntStarts - A Float64 vector that gives the starting points of intervals (in the x plane)
  * SpCoefs - A 3 column matrix with the same number of rows as the length of the IntStarts vector. The first column is the coefficient of the quadratic term, the second column for the linear term. The third column is the constant.

**Returns**

  * A spline function that takes a single Float64 input and returns the spline value at that point.  

<a id='SchumakerSpline.ppmakDeriv' href='#SchumakerSpline.ppmakDeriv'>#</a>
**`SchumakerSpline.ppmakDeriv`** &mdash; *Function*.



Creates the derivative function of the spline defined by interval starts IntStarts and quadratic coefficients SpCoefs which evaluates an input point.

**Takes**

  * IntStarts - A Float64 vector that gives the starting points of intervals (in the x plane)
  * SpCoefs - A 3 column matrix with the same number of rows as the length of the IntStarts vector. The first column is the coefficient of the quadratic term, the second column for the linear term. The third column is the constant.

**Returns**

  * The derivative function that takes a single Float64 input and returns the derivative at that point.  

<a id='SchumakerSpline.ppmak2Deriv' href='#SchumakerSpline.ppmak2Deriv'>#</a>
**`SchumakerSpline.ppmak2Deriv`** &mdash; *Function*.



Creates the second derivative function of the spline defined by interval starts IntStarts and quadratic coefficients SpCoefs which evaluates an input point.

**Takes**

  * IntStarts - A Float64 vector that gives the starting points of intervals (in the x plane)
  * SpCoefs - A 3 column matrix with the same number of rows as the length of the IntStarts vector. The first column is the coefficient of the quadratic term, the second column for the linear term. The third column is the constant.

**Returns**

  * The second derivative function that takes a single Float64 input and returns the second derivative at that point.

