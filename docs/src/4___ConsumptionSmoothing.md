# 4 Applications
## 4.4 A consumption smoothing problem

Consider an infinitely lived consumer that has a budget of $B_t$ at time $t$ and a periodic income of $1$. She has a periodic utility function given by $\epsilon_t x_t^\delta$, where $x_t$ is spending in period $t$ and $\epsilon_t$ is the shock in period $t$ drawn from some stationary nonnegative shock process with pdf $f(\epsilon)$ defined on the interval $[y,z]$. The problem for the consumer in period $t$ is to maximise their value function:

$V(B_t | \epsilon_{t}) =  \max_{0 < x_t < B_t} \hspace{0.5cm} \epsilon_t x_t^\delta + \beta \int_y^z V(B_{t+1} | \epsilon) f(\epsilon)d\epsilon$

Where $\beta$ is a discounting factor and $B_{t+1} = 1 + B_t - x_t$.

Our goal is that we want to find a function that gives the optimal spending amount, $\hat{x}(B_t, \epsilon_t)$,  in period $t$ which is a function of the shock magnitude $\epsilon_{t}$ and the saved budgets $B_{t}$ in this period. If we knew the function $\int_y^z V(B_{t+1} \vert \epsilon) f(\epsilon)d\epsilon$ then we could do this by remembering $B_{t+1} = 1 + B_t - x_t$ and using the optimisation:

$\hat{x}(B_t, \epsilon_t) = \text{argmax}_{0 < x_t < B_t} \hspace{0.5cm} \epsilon_t x_t^\delta + \beta \int_y^z V(B_{t+1} | \epsilon) f(\epsilon)d\epsilon$

So now we need to find the function $E_t[ V(B_{t+1})]$. Note as the shock process is stationary, the consumer lives forever and income is always 1, $E_t[ V(B_{t+1})]$ will not vary with $t$. As a result we will rewrite it as simply $f(b)$.

Now we will construct a vector containing a grid of budget values, $\bar{b}$, for instance $\bar{b} = [0, 0.01,0.02, ... , 5]$ (we will use bars to describe approximations gained from this grid). If we could then approximate a vector of the corresponding function values, $\bar{f}$,  so we had for instance $\bar{f} = [f(0), f(0.01), f(0.02), ... , f(5)]$ then we could approximate the function by constructing a spline $\bar{f}(b)$ between these points. Then we can get the function:

$\bar{x}(B_t, \epsilon_t) = \text{argmax}_{0 < x < B_t} \hspace{0.5cm} \epsilon_t x_t^{\delta} + \bar{f}(B_{t} - x)]$

So this problem reduces to finding the vector of function values at a discrete number of points, $\bar{f}$. This can be done as a fixed point problem. We can first note that this problem is a contraction mapping problem. In this particular example this means that if we define a sequence $\bar{f}_0 = f_0$ where $f_0$ is some initial guess and $f_i = g(f_{i-1})$ where $g$ is given by the IterateOnce function below then this sequence will be convergent. Convergence would be slow however so below we will actually use the Anderson method:


```
using Distributions
using FixedPointAcceleration
using HCubature
using Optim
using Random
using SchumakerSpline
delta = 0.2
beta = 0.95
periodic_income = 1.0
shock_var = 1.0
shock_process = LogNormal(0.0, shock_var)
BudgetStateSpace = vcat( collect(0:0.015:periodic_income), collect(1.05:0.05:(3*periodic_income)))
InitialGuess = sqrt.(BudgetStateSpace)

function ValueGivenShock(Budget::Float64, epsilon::Float64, NextValueFunction::Schumaker)
    opt = optimize(x ->  -1.0*(epsilon*(x^delta) + beta*evaluate(NextValueFunction, Budget - x + periodic_income)), 0.0, Budget)
    return -1.0 * opt.minimum
end

function ExpectedUtility(Budget::Float64, NextValueFunction::Schumaker)
    if Budget > 0.00001
        integ = hcubature(epsilon-> ValueGivenShock(Budget, epsilon[1], NextValueFunction)* pdf(shock_process, epsilon[1]), [quantile(shock_process,0.0001)], [quantile(shock_process, 0.9999)])
        return integ[1]
    else
        return beta * evaluate(NextValueFunction, periodic_income)
    end
end

function OneIterateBudgetValues(BudgetValues::Array{Float64,1})
    NextValueFunction = Schumaker(BudgetStateSpace, BudgetValues)
    new_budget_values = zeros(length(BudgetStateSpace))
    for i in 1:length(BudgetStateSpace)
        new_budget_values[i] = ExpectedUtility(BudgetStateSpace[i], NextValueFunction)
    end
    return new_budget_values
end

fp_anderson = fixed_point(OneIterateBudgetValues, InitialGuess; Algorithm = Anderson, PrintReports = true)
fp_simple   = fixed_point(OneIterateBudgetValues, InitialGuess; Algorithm = Simple, PrintReports = true)
```

This takes 22 iterates with the Anderson algorithm which is drastically better than the 459 iterates it takes with the simple method.
