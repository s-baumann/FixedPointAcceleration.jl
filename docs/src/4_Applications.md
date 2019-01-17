# 4.0 Applications

## 4.1 Finding equilibrium prices in a pure exchange economy

Consider that there are N households in a pure exchange economy. Every household has preferences over G types of good. Household $n$ has a utility function of

$U_n = \sum_{i=1}^G \gamma_{n,i} \log(c_{n,i})$

Where $\gamma_{n,i}$ is a parameter describing household $n$'s taste for good $i$, $c_{n,i}$ is
household $n$'s consumption of good $i$.
Each household is endowed with an amount of each good. They can then trade goods before consumption.
We want to find the equilibrium prices in this exchange economy. We have data on each household's
endowment and preferences for each good and want to determine the equilibrium prices for this pure
exchange economy.

We will choose good 1 as the numeraire. So $P_1 = 1$. First we will find an expression for demand given a price vector. Setting up the lagrangian for household $n$:

$L_n = \sum_{i=1}^G \gamma_{n,i} \log(c_{n,i}) + \lambda_{n}[ \sum_{i=1}^G p_i(e_{n,i}-c_{n,i}) ]$

Where $\lambda_{n}$ is household $n$'s shadow price and $e_{n,i}$ is this household's endowment of
good $i$ and $p_i$ is the price of good $i$. Taking FOC with respect to $c_i$ of this lagrangian yields:

$c_{n,i} = \frac{\gamma_{n,i}}{p_i \lambda_n}$

and taking FOC condition with respect to $\lambda_n$ yields the budget constraint. Subbing the above
equation into the budget constrain and rearranging yields:

$\lambda_n = \frac{\sum^G_{i=1} \gamma_{n,i}}{\sum^G_{i=1} p_i e_{n,i}}$

We can also sum over households to find total demand for each good as:

$D_i = \frac{1}{P_i} \sum_{n=1}^G \frac{\gamma_{n,i}}{\lambda_n}$

We will find the equilibrium price vector by using an approximate price vector to find the $\lambda$s. We
can then find an estimate of the equilibrium price $P_i$ which solves $D_i = \sum_{n=1}^G e_{n,i}$:

$P_i = \frac{\sum_{n=1}^G e_{n,i}}{\sum_{n=1}^G \frac{\gamma_{n,i}}{\lambda_n} }$

We use this approach in the code below for the case of 10 goods with 8 households. For exposition sake we generate some data below before proceeding to find the equilibrium price vector.

```
# Generating data
using Distributions
using FixedPointAcceleration
using Random
Random.seed!(1234)
N = 5
G = 10
Endowments = rand(LogNormal(), G, N)
Tastes      = rand(G, N)  
# Every column here represents a household and every row is a good. So Endowments[1,2] is
# the second household's endowment of good 1.

# We now start solving for equilibrium prices:
TotalEndowmentsPerGood = mapslices(sum, Endowments; dims = [2])
TotalTastesPerHousehold = mapslices(sum, Tastes; dims = [1])

function LambdasGivenPriceVector(prices)
    ValueOfEndowmentsPerHousehold = prices .* Endowments
    TotalValueOfEndowmentsPerHousehold =  mapslices(sum, ValueOfEndowmentsPerHousehold; dims = [1])
    return TotalTastesPerHousehold ./ TotalValueOfEndowmentsPerHousehold
end

function IterateOnce(prices)
    Lambdas = LambdasGivenPriceVector(prices)
    TastesOverLambdas = Tastes ./ Lambdas
    SumTastesOverLambdas = mapslices(sum, TastesOverLambdas; dims = [2])
    NewPrices = SumTastesOverLambdas ./ TotalEndowmentsPerGood
    NewPrices = NewPrices/NewPrices[1] # Applying Numeraire
    return NewPrices
end


InitialGuess = repeat([1.0], 10)
FPSolution = fixed_point(IterateOnce, InitialGuess; Algorithm = VEA)
```

## 4.2 A consumption smoothing problem

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

This takes 22 iterates with the anderson algorithm which is drastically better than the several hundred iterates it takes with the simple method.