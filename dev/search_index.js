var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "FixedPointAcceleration.jl",
    "title": "FixedPointAcceleration.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#FixedPointAcceleration.jl-1",
    "page": "FixedPointAcceleration.jl",
    "title": "FixedPointAcceleration.jl",
    "category": "section",
    "text": "Fixed point finders are conceptually similar to both optimisation and root finding algorithms but thus far implementations of fixed point finders have been rare in Julia. In some part this is likely because there is often an obvious method to find a fixed point by merely feeding a guessed fixed point into a function, taking the result and feeding it back into the function. By doing this repeatedly a fixed point is often found. This method (that we will call the \"Simple\" method) is often convergent but it is also often slow which can be prohibitive when the function itself is expensive.FixedPointAcceleration.jl aims to provide fixed point acceleration algorithms that can be much faster than the simple method. In total eight algorithms are implemented. The first is the simple method as described earlier. There are also the Newton, Aitken and Scalar Epsilon Algorithm (SEA) methods that are designed for accelerating the convergence of scalar sequences. Four algorithms for accelerating vector sequences are also implemented including the Vector Epsilon Algorithm (VEA), two minimal polynomial algorithms (MPE and RRE)  and Anderson acceleration.In this paper section 1 starts by with a brief explanation of fixed points before section 2 describes the acceleration algorithms provided by FixedPointAcceleration.jl. Here the goal is  to illustrate the logic underling each algorithm so users can better choose which suits their problem. Readers interested in the underlying mathematics are referred to the original papers. Section 3 then illustrates a few features of the package that enable users to better track the progress of an algorithm while it is running and switch algorithms if desired before a fixed point is found.Section 4 then presents several applications of these fixed point algorithms in economics, asset pricing and machine learning. Finally section 5 presents a convergence comparison showing how many iterations each algorithm takes in bringing each problem to its fixed point for each of the applications presented in section 4."
},

{
    "location": "#Fixed-point-acceleration-1",
    "page": "FixedPointAcceleration.jl",
    "title": "1 Fixed point acceleration",
    "category": "section",
    "text": "A fixed point problem is one where we look for a vector, hatX in Re^N, so that for a given function f Re^N rightarrow Re^N we have:f(hatX) = hatXIf f Re^1 rightarrow Re^1 and thus any solution hatX will be a scalar then one way to solve this problem would be to use a rootfinder on the function g(x) = f(x) - x or to use an optimiser to minimise h(x) = (f(x) - x)^2. These techniques will not generally work however if f  N^a rightarrow N^a where a is large. Consider for instance using a multidimensional Newtonian optimiser to minimise h(x) = (f(x) - x)^2. The estimation of gradients for each individual dimension may take an infeasibly long time. In addition this method may not make use all available information. Consider for instance that we know that the solution for x will be an increasing vector (so x_i  x_j for any entries of x  with i  j) but has many entries. This information can be preserved and used in the vector acceleration algorithms that we present but would be more difficult to exploit in a standard optimisation algorithm.FixedPointAcceleration.jl implements eight algorithms for finding fixed points. The first algorithm implemented in this package is the \"simple\" method which merely takes the output of a function and feeds it back into the function. For instance starting with a guess of x_0, the next guess will be x_1 = f(x_0). The guess after that will be x_2 = f(x_1) and so on. In some conditions f will be a contraction mapping and so the simple method will be guaranteed to converge to a unique fixed point (Stokey, Lucas & Prescott 1989). Even when this is the case however the simple method may only converge slowly which can be inconvenient. The other seven methods this package implements are designed to be faster than the simple method but may not be convergent for every problem.pages = [\"index.md\",\n         \"2_Algorithms.md\",\n         \"3_UsingAdvice.md\",\n         \"99_refs.md\"]\nDepth = 2"
},

{
    "location": "2_Algorithms/#",
    "page": "2 Acceleration algorithms",
    "title": "2 Acceleration algorithms",
    "category": "page",
    "text": ""
},

{
    "location": "2_Algorithms/#Acceleration-algorithms-1",
    "page": "2 Acceleration algorithms",
    "title": "2 Acceleration algorithms",
    "category": "section",
    "text": ""
},

{
    "location": "2_Algorithms/#.1-Newton-acceleration-1",
    "page": "2 Acceleration algorithms",
    "title": "2.1 Newton acceleration",
    "category": "section",
    "text": "Here we will define g(x) = f(x) - x. The general approach is to solve g(x) with a rootfinder. The x that provides this root will be a fixed point. Thus after two iterates we can approximate the fixed point with:x_i+1 = x_i - fracg(x_i)g^prime(x_i)FixedPointAcceleration.jl approximates the derivative g^prime(x_i) so that we use an estimated fixed point of:x_i+1 = x_i - fracg(x_i)  frac g(x_i) - g(x_i-1)x_i-x_i-1      The implementation of the Newton method in this package uses this formula to predict the fixed point given two previous iterates.[1] This method is designed for use with scalar functions.[1]: Only two rather than three are required because we can use x_i+1 = f(x_i) and x_i+2 = f(x_i+1)."
},

{
    "location": "2_Algorithms/#.2-Aitken-acceleration-1",
    "page": "2 Acceleration algorithms",
    "title": "2.2 Aitken acceleration",
    "category": "section",
    "text": "Consider that a sequence of scalars  x_i _i=0^infty that converges linearly to its fixed point of hatx. This implies that for a large i:frachatx - x_i+1hatx - x_i approx frachatx - x_i+2hatx - x_i+1For a concrete example consider that every iteration halves the distance between the current vector and the fixed point. In this case the left hand side will be one half which will equal the right hand side which will also be one half.This above expression can be simply rearranged to give a formula predicting the fixed point that is used as the subsequent iterate. This is:textNext Guess = x_i - frac  (x_i+1 - x_i)^2    x_i+2 - 2x_i+1 + x_iThe implementation of the Aitken method in FixedPointAcceleration.jl uses this formula to predict the fixed point given two previous iterates. This method is designed for use with scalar functions. If it is used with higher dimensional functions that take and return vectors then it will be used elementwise."
},

{
    "location": "2_Algorithms/#.3-Epsilon-algorithms-1",
    "page": "2 Acceleration algorithms",
    "title": "2.3 Epsilon algorithms",
    "category": "section",
    "text": "The epsilon algorithms introduced by Wynn (1962) provides an alternate method to extrapolate to a fixed point. This documentation will present a brief numerical example and refer readers to Wynn (1962) or  Smith, Ford & Sidi (1987) for a mathematical explanation of why it works. The basic epsilon algorithm starts with a column of iterates (column B in the below figure). If i iterates have been performed then this column will have a length of i+1 (the initial starting guess and the results of the i iterations). Then a series of columns are generated by means of the below equation:epsilon_c+1^r = epsilon_c-1^r+1 + (epsilon_c^r+1 - epsilon_c^r)^-1Where c is a column index and r is a row index. The algorithm starts with the epsilon_0 column being all zeros and epsilon_1 being the column of the sequence iterates. The value in the furthest right column ends up being the extrapolated value.This can be seen in the below table which uses an epsilon method to find the fixed point of cos(x) with an initial guess of a fixed point of 1.(Image: The Epsilon Algorithm applied to the cos(x) function)In this figure B1 is the initial guess of the fixed point. Then we have the iterates B2 = cos(B1), B3 = cos(B2) and so on. Moving to the next column we have C1 = A2 + 1(B2-B1) and C2 = A3 + 1(B3-B2) and so on before finally we get F1 = D2 + 1(E2-E1). As this is the last entry in the triangle it is also the extrapolated value.Note that the values in columns C and E are poor extrapolations. Only the even columns D,F provide reasonable extrapolation values. For this reason an even number of iterates (an odd number of values including the starting guess) should be used for extrapolation. FixedPointAcceleration.jl\'s functions will enforce this by throwing away the first iterate provided if necessary to get an even number of iterates.In the vector case this algorithm can be visualised by considering each entry in the above table to contain a vector going into the page. In this case the complication emerges from the inverse - there is no clear interpretation of (epsilon_c^r+1 - epsilon_c^r)^-1 where (epsilon_c^r+1 - epsilon_c^r) represents a vector. The Scalar Epsilon Algorithm (SEA) uses elementwise inverses to solve this problem which ignores the vectorised nature of the function. The Vector Epsilon Algorithm (VEA) uses the Samuelson inverse of each vector (epsilon_c^r+1 - epsilon_c^r)."
},

{
    "location": "2_Algorithms/#.4-Minimal-polynomial-algorithms-1",
    "page": "2 Acceleration algorithms",
    "title": "2.4 Minimal polynomial algorithms",
    "category": "section",
    "text": "FixedPointAcceleration.jl implements two minimal polynomial algorithms, Minimal Polynomial Extrapolation (MPE) and Reduced Rank Extrapolation (RRE). This section will simply present the main equations but an interested reader is directed to Cabay & Jackson (1976) or Smith, Ford & Sidi (1987) for a detailed explanation.To first define notation, each vector (the initial guess and subsequent iterates) is defined by x_0 x_1  . The first differences are denoted u_j = x_j+1 - x_j and the second differences are denoted v_j = u_j+1 - u_j. If we have already completed k-1 iterations (and so we have k terms) then we will use matrices of first and second differences with U =  u_0  u_1    u_k-1  and V =  v_0  v_1   v_k-1 .For the MPE method the extrapolated vector, s, is found by:s = frac sum^k_ j=0   c_j x_j  sum^k_j=0 c_j Where the coefficient vector is found by c = -U^+ u_k where U^+  is the Moore-Penrose generalised inverse of the U matrix. In the case of the RRE the extrapolated vector, s, is found by:s = x_0 - U V^+ u_0"
},

{
    "location": "2_Algorithms/#.5-Anderson-acceleration-1",
    "page": "2 Acceleration algorithms",
    "title": "2.5 Anderson acceleration",
    "category": "section",
    "text": "Anderson (1965) acceleration is an acceleration algorithm that is well suited to functions of vectors. Similarly to the minimal polynomial algorithms it takes a weighted average of previous iterates. It is different however to these algorithms (and the VEA algorithm) in that previous iterates need not be sequential but any previous iterates can be used.Consider that we have previously run an N-dimensional function M times. We can define a matrix G_i = g_i-M  g_i-M+1   g_i where g(x_j) = f(x_j) - x_j. Each column of this matrix can be interpreted as giving the amount of ``movement\'\' that occurred in a run of the function.Now Anderson acceleration defines a weight to apply to each column of the matrix. This weight vector is M-dimensional and can be denoted mathbfalpha = alpha_0 alpha_1    alpha_M. These weights are determined by means of the following optimisation:min_mathbfalpha vertvert G_i mathbfalpha vertvert_2hspace1cm st sum^M_j=0 alpha_j = 1Thus we choose the weights that will be predicted to create the lowest ``movement\'\' in an iteration.With these weights we can then create the expression for the next iterate as:x_i+1 = sum_j=0^M alpha_j f(x_i-M+j)The actual implementation differs from the proceeding description by recasting the optimisation problem as an unconstrained least squares problem (see Fang & Saad 2009 or Walker & Ni 2011) but in practical terms is identical."
},

{
    "location": "3_UsingAdvice/#",
    "page": "3.  Using the FixedPoint package",
    "title": "3.  Using the FixedPoint package",
    "category": "page",
    "text": ""
},

{
    "location": "3_UsingAdvice/#.-Using-the-FixedPoint-package-1",
    "page": "3.  Using the FixedPoint package",
    "title": "3.  Using the FixedPoint package",
    "category": "section",
    "text": ""
},

{
    "location": "3_UsingAdvice/#.1-Basic-examples-of-using-FixedPoint-1",
    "page": "3.  Using the FixedPoint package",
    "title": "3.1 Basic examples of using FixedPoint",
    "category": "section",
    "text": ""
},

{
    "location": "3_UsingAdvice/#The-Babylonian-method-for-finding-square-roots.-1",
    "page": "3.  Using the FixedPoint package",
    "title": "The Babylonian method for finding square roots.",
    "category": "section",
    "text": "Now we will demonstrate how FixedPoint can be used for simple problems. For the simplest possible case consider we want to estimate a square root using the Babylonian method. To find the square root of a number x, given an initial guess t_0 the following sequence converges to the square root:t_n+1 = frac12 left t_n + fracxt_n rightThis is a fast converging and inexpensive sequence which probably makes an acceleration algorithm overkill but for sake of exposition we can implement this in FixedPoint. In the next code block we find the square root of 100 with the simple method:using FixedPointAcceleration\nSequenceFunction = function(tn){0.5*(tn + 100/tn)}\nInitial_Guess = 6.0\nFP_Simple   = fixed_point(SequenceFunction,Initial_Guess; Method = \"Simple\")We can also solve for a vector of fixed points at the same time. For instance every square root from 1 to 100.NumbersVector = 1:100\nSequenceFunction = function(tn){0.5*(tn + NumbersVector/tn)}\nFP_SEA   = FixedPoint(Function = SequenceFunction, Inputs = 1:100,  Method = \"RRE\")Note that in this case the RRE method is being applied elementwise."
},

{
    "location": "3_UsingAdvice/#Vectorised-functions-1",
    "page": "3.  Using the FixedPoint package",
    "title": "Vectorised functions",
    "category": "section",
    "text": "The utility of the acceleration algorithms contained in FixedPoint are more apparent when applied to vectorised functions with cross dependency. For a simple example consider the below function where each entry of the vector depends on both entries of the previous iterate.SimpleVectorFunction = function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])}\nFP_Simple = FixedPoint(Function = SimpleVectorFunction, Inputs = c(0.3,900), Method = \"Simple\")\nFP_Anderson = FixedPoint(Function = SimpleVectorFunction, Inputs = c(0.3,900), Method = \"Anderson\")This function takes 105 iterates to find a fixed point with the simple method but only 14 with the Anderson acceleration method."
},

{
    "location": "3_UsingAdvice/#.2-Easily-changing-algorithm-1",
    "page": "3.  Using the FixedPoint package",
    "title": "3.2 Easily changing algorithm",
    "category": "section",
    "text": "A key priority in writing FixedPoint is for ease in changing fixed point method. The algorithm can be stopped and then restarted without having to repeat earlier iterations. This can be seen in the simple example below:Run1 = FixedPoint(Function = function(x){x^0.2}, Inputs = c(0.3, 0.4,0.5),\n                Method = \"Aitken\")\ncat(paste(\"Running this all at once takes \", length(Run1$Convergence), \" iterations and\n          convergences to \",  Run1$Convergence[length(Run1$Convergence)], \"\\n\"))## Running this all at once takes  9  iterations and\n##           convergences to  1.11022302462516e-15Run2_1 = FixedPoint(Function = function(x){x^0.2}, Inputs = c(0.3, 0.4,0.5),\n                Method = \"Aitken\", MaxIter = 4)\n# Check some condition about the convergence or change algorithm.\nRun2_2 = FixedPoint(Function = function(x){x^0.2}, Inputs = Run2_1$Inputs,\n                    Outputs = Run2_1$Outputs, Method = \"Aitken\")\ncat(paste(\"Running this with a break takes \", length(Run2_2$Convergence), \" iterations and\n          convergences to \",  Run2_2$Convergence[length(Run2_2$Convergence)], \"\\n\"))## Running this with a break takes  9  iterations and\n##           convergences to  1.11022302462516e-15This can be useful in problems where it is desired to switch extrapolation methods half way through. For instance in the consumption smoothing case presented later, the Vector Epsilon Algorithm converges quickly to a convergence level around 10^-4 but then spends around 100 iterations hovering around this level without reducing further. In this case the algorithm could be changed to the simple method or Anderson method once the algorithm hit some target level of convergence.Another example is parallelisation. The Anderson method in particular is well suited to use the results of previous iterates and these iterates do not need to be generated sequentially. If a suitable collection of guess vectors can be generated[4] then all could be evaluated in parallel and recombined before again being used by the Anderson algorithm to generate a new guess. An example of this is presented later when the package is used to find the value function for a consumer in a consumption smoothing problem.[4]: For instance a set of guesses could be done by using different dampening values or by running the Anderson algorithm several times using different subsets of previous iterates to generate new guesses."
},

{
    "location": "3_UsingAdvice/#.4-Graceful-error-handling-1",
    "page": "3.  Using the FixedPoint package",
    "title": "3.4 Graceful error handling",
    "category": "section",
    "text": "Hopefully FixedPoint is well tested enough that most kind of errors will be rare. FixedPoint also offers an option (ReplaceInvalids) to ensure that no acceleration algorithm generates guess vectors that contain NAs or Infs This option can be set toReplaceVector\'\'  which will replace an extrapolated vector containing NAs or Infs by the vector output in the previous iterate. If it is set to ``ReplaceElement\'\' then it will replace the individual elements that are NA or Inf by the corresponding elements in the output of the previous iterate.[15][15]: With these options being set the FixedPoint should never throw an error when applied to a function which is guaranteed to map: Re^N rightarrow Re^N for some N unless the error comes from underflow or overflow. If you do encounter such a testcase it likely reflects a bug, we would much appreciate it if you could send a minimal working example to the package maintainer.Errors are likely however in cases where inputs functions have a restricted domain. For example this may include functions that require the input vector to have a particular shape (ie concavity) or functions where the input vector must be strictly positive. For a simple example consider the vectorised function we introduced in section 3.1. Now rather than x^prime1 = fracsqrtvert x1 + x2 vert2 we have x^prime1 = fracsqrt x1 + x2 2 where the output x has a prime and the inputs has no prime symbol. x^prime1 here is no longer real valued if the sum of the previous iterate is negative. This is what occurs in the 5th iterate of the Anderson method applied to this problem.The FixedPoint function handles these situations gracefully by saving all previous results as well as the proposed new vector that lead to the error. In these cases the \"Finish\" message in the output describes the form of error, the \"NewInputVector\" list entry contains the attempted vector that was input to the function and the \"NewOutputVector\" list entry contains the result of applying the function to NewInputVector. This output information can be seen in the example below.This information is useful in order to diagnose the issue. In this case we might decide to modify the function to insert the absolute value function with the reasoning that the same fixed point will likely result (which we could later verify). This also allows a user to run one method until an error occurs and then switch methods. This is demonstrated below.SimpleVectorFunction = function(x){c(0.5*sqrt(x[1] + x[2]), 1.5*x[1] + 0.5*x[2])}\nFPSolution = FixedPoint(Function = SimpleVectorFunction, Inputs = c(0.3,900),\n                        Method = \"Anderson\")## Warning in sqrt(x[1] + x[2]): NaNs produced# We can see the output of the FixedPoint function in cases like this where it ends due\n# to an error\nFPSolution## $Inputs\n##       [,1]     [,2]      [,3]     [,4]\n## [1,]   0.3  15.0025  7.350817 3.191655\n## [2,] 900.0 450.4500 82.469248 9.574966\n##\n## $Outputs\n##          [,1]      [,2]      [,3]     [,4]\n## [1,]  15.0025  10.78717  4.738672 1.786520\n## [2,] 450.4500 247.72875 52.260849 9.574966\n##\n## $Convergence\n## [1] 449.550000 202.721250  30.208399   1.405135\n##\n## $FixedPoint\n## [1] NA\n##\n## $Finish\n## [1] \"New output vector contains NAs\"\n##\n## $NewInputVector\n## [1] -1.085113 -3.255338\n##\n## $NewOutputVector\n## [1]       NaN -3.255338# We can use this information to decide to switch to the simple method.\n# No error results as the simple method doesn\'t extrapolate.\nFPSolution = FixedPoint(Function = SimpleVectorFunction, Inputs = FPSolution$Inputs,\n                        Outputs = FPSolution$Outputs, Method = \"Simple\", MaxIter = 5)\n# Now we switch to the Anderson Method again. No error results because we are\n# close to fixed point.\nFPSolution = FixedPoint(Function = SimpleVectorFunction, Inputs = FPSolution$Inputs,\n                        Outputs = FPSolution$Outputs, Method = \"Anderson\")"
},

{
    "location": "3_UsingAdvice/#.5-Convergence-by-constant-increments-1",
    "page": "3.  Using the FixedPoint package",
    "title": "3.5 Convergence by constant increments",
    "category": "section",
    "text": "Most of the methods included in this function will fail in finding the fixed point of a function that converges by a fixed increment. For instance we may have a function that takes x and returns x shifted 1 unit (in Euclidian norm) in a straight line towards its fixed point. A realistic example of this is the training of a perceptron classifier which is explored later in section 4.3.This case is problematic for all methods except for the simple method. The basic problem can be illustrated simply by looking at the Newton method and the Aitken method. For the Newton method the derivative is approximated by frac g(x_i) - g(x_i-1)x_i-xi-1. When there is convergence by constant increments then g(x_i) = g(x_i-1)  and the derivative is zero which means calculating the Newton methods recommended new guess of the fixed point involves division by zero Now considering the Aitken method the new guess is given by x_i+1 = x_i - frac  (x_i+1 - x_i)^2    x_i+2 - 2x_i+1 + x_i. When there is convergence by constant increments then x_i - x_i+1 = x_i+1 - x_i+2  and so we have x_i+2 - 2x_i+1 + x_i = (x_i - x_i+1) - (x_i+1 - x_i+2) = 0. It is against not possible to calculate the new guess.[5]More generally similar problems exist for the other acceleration methods. When there is convergence by constant increments then then the fixed point method receives information about what direction to go in but no information about how far to go. This is a complication that is common to all of these acceleration methods in this package. In these cases it may be possible to change the function to make it converge by varying increments while retaining the same set of fixed points. This is shown in the perceptron example in section 4.3. In other cases where it is not possible to modify the function, it is advisable to use the simple method.[5]: When these complications arise the ReplaceInvalids method can be used to revert to a simple iterate or to change individual elements to the corresponding values in a simple iterate. This is as described in section 3.4."
},

{
    "location": "99_refs/#",
    "page": "References",
    "title": "References",
    "category": "page",
    "text": ""
},

{
    "location": "99_refs/#References-1",
    "page": "References",
    "title": "References",
    "category": "section",
    "text": "Anderson, D.G. 1965. “Iterative Procedures for Nonlinear Integral Equations.” Journal of the ACM 12 (4): 547–60.Cabay, S., and L.W. Jackson. 1976. “A Polynomial Extrapolation Method for Finding Limits and Antilimits of Vector Sequences.” Siam Journal of Numerical Analysis 13 (5): 734–52.Fang, Haw-ren, and Yousef Saad. 2009. “Two Classes of Multisecant Methods for Nonlinear Acceleration.” Numerical Linear Algebra with Applications 16 (3): 197–221.Smith, David, William Ford, and Avram Sidi. 1987. “Extrapolation Methods for Vector Sequences.” SIAM Review 29 (2): 199–233.Stokey, Nancy, Robert E. Lucas, and Edward Prescott. 1989. Recursive Methods in Economic Dynamics. Harvard University Press.Walker, Homer. 2010. “Anderson Acceleration: Algorithms and Implementations.” https://users.wpi.edu/~walker/Papers/andersonaccnalgs_imps.pdf.Walker, Homer, and Peng Ni. 2011. “Anderson Acceleration for Fixed-Point Iterations.” SIAM Review 49(4): 1715–35.Wynn, P. 1962. “Acceleration Techniques for Iterated Vector and Matrix Problems.” Mathematics of Computation 16 (79): 301–22."
},

]}
