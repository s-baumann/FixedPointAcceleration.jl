# FixedPointAcceleration.jl

Fixed point finders are conceptually similar to both optimisation and root finding algorithms but thus far implementations of fixed point finders have been rare in Julia. In some part this is likely because there is often an obvious method to find a fixed point by merely feeding a guessed fixed point into a function, taking the result and feeding it back into the function. By doing this repeatedly a fixed point is often found. This method (that we will call the "Simple" method) is often convergent but it is also often slow which can be prohibitive when the function itself is expensive.

**FixedPointAcceleration.jl** aims to provide fixed point acceleration algorithms that can be much faster than the simple method. In total eight algorithms are implemented. The first is the simple method as described earlier. There are also the Newton, Aitken and Scalar Epsilon Algorithm (SEA) methods that are designed for accelerating the convergence of scalar sequences. Four algorithms for accelerating vector sequences are also implemented including the Vector Epsilon Algorithm (VEA), two minimal polynomial algorithms (MPE and RRE)  and Anderson acceleration.

In this paper section 1 starts by with a brief explanation of fixed points before section 2 describes the acceleration algorithms provided by **FixedPointAcceleration.jl**. Here the goal is  to illustrate the logic underling each algorithm so users can better choose which suits their problem. Readers interested in the underlying mathematics are referred to the original papers. Section 3 then illustrates a few features of the package that enable users to better track the progress of an algorithm while it is running and switch algorithms if desired before a fixed point is found.

Section 4 then presents several applications of these fixed point algorithms in economics, asset pricing and machine learning. Finally section 5 presents a convergence comparison showing how many iterations each algorithm takes in bringing each problem to its fixed point for each of the applications presented in section 4.

# 1 Fixed point acceleration

A fixed point problem is one where we look for a vector, $\hat{X} \in \Re^N$, so that for a given function $f: \Re^N \rightarrow \Re^N$ we have:

$$f(\hat{X}) = \hat{X}$$

If $f: \Re^1 \rightarrow \Re^1$ and thus any solution $\hat{X}$ will be a scalar then one way to solve this problem would be to use a rootfinder on the function $g(x) = f(x) - x$ or to use an optimiser to minimise $h(x) = (f(x) - x)^2$. These techniques will not generally work however if $f : N^a \rightarrow N^a$ where $a$ is large. Consider for instance using a multidimensional Newtonian optimiser to minimise $h(x) = (f(x) - x)^2$. The estimation of gradients for each individual dimension may take an infeasibly long time. In addition this method may not make use all available information. Consider for instance that we know that the solution for $x$ will be an increasing vector (so $x_i > x_j$ for any entries of $x$  with $i > j$) but has many entries. This information can be preserved and used in the vector acceleration algorithms that we present but would be more difficult to exploit in a standard optimisation algorithm.

**FixedPointAcceleration.jl** implements eight algorithms for finding fixed points. The first algorithm implemented in this package is the "simple" method which merely takes the output of a function and feeds it back into the function. For instance starting with a guess of $x_0$, the next guess will be $x_1 = f(x_0)$. The guess after that will be $x_2 = f(x_1)$ and so on. In some conditions $f$ will be a contraction mapping and so the simple method will be guaranteed to converge to a unique fixed point (Stokey, Lucas & Prescott 1989). Even when this is the case however the simple method may only converge slowly which can be inconvenient. The other seven methods this package implements are designed to be faster than the simple method but may not be convergent for every problem.

```@contents
pages = ["index.md",
         "2_Algorithms.md",
         "3_UsingAdvice.md",
         "99_refs.md"]
Depth = 2
```
