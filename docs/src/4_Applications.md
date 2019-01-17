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

## 4.2 The Perceptron Classifier

The perceptron is one of the oldest and simplest machine learning algorithms (Rosenblatt 1958). In its simplest form, for each observation it is applied it uses an N-dimensional vector of features x together with N+1 weights w to classify the observation as being in category one or category zero. It classifies observation j as a type one if $$w_0 + \sum_{i=1}^N w_i x_{i,j}  > 0$$ and as a type zero otherwise.

The innovation of the perceptron was its method for training its weights, w. This is done by looping over a set of observations that can be used for training (the "training set") and for which the true category information is available.
The perceptron classifies each observation. When it correctly classifies an observation no action is taken. On the other hand when the perceptron makes an error then it updates its weights with the following expressions.

$$w_{0}^\prime = w_{0} + ( d_{j} - y_{j} )$$

$$w_{i}^\prime = w_{i} + ( d_{j} - y_{j} ) x_{j,i} \hspace{1cm} \text{ for } i \geq 0$$

Where $$w_i$$ is the old weight for the $i$'th feature and $$w_{i}^\prime$$ is the updated weight. $$x_{j,i}$$ is the feature value for
observation $j$'s feature $i$, $$d_{j}$$ is the category label for observation $j$ and $$y_j$$ is the perceptron's prediction for this
observation’s category.

This training algorithm can be rewritten as fixed point problem. We can write a function that takes perceptron weights, loops over the
data updating these weights and then returns the updated weight vector. If the perceptron classifies every observation correctly then
the weights will not update and we are at a fixed point.[^7]

The standard Perceptron updating algorithm does not work well with FixedPointAcceleration methods because of convergence by a fixed increment. This occurs because multiple iterates can result in the same observations being misclassified and hence the same change in the weights. As a result we modify the training algorithm to give training increments that change depending on distance from the fixedpoint. This can be done by updating the weights by an amount proportional to a concave function of the norm of $wx+b$.

[^7]: Note that for perceptrons there are always uncountably many such fixed points
where the perceptron correctly classifies the entire training set and will not further update. On the other hand it is possible that
the data is not linearly separable in which case there may be no fixed point and the weights will continue to update forever.

First we generate a dataset:
```
# Generating linearly seperable data
using Distributions
using FixedPointAcceleration
using Random
using DataFrames
nobs = 20
Random.seed!(1234)
data1 = DataFrame([rand(Normal(3,2), nobs), rand(Normal(8,2), nobs), repeat([-1.0],nobs)], [:x1, :x2, :y])
data2 = DataFrame([rand(Normal(-4,2), nobs), rand(Normal(10,12), nobs), repeat([1.0],nobs)], [:x1, :x2, :y])
data  = vcat(data1,data2)
# Plotting it
using Plots
plot(data1.x1, data1.x2,seriestype=:scatter)
plot!(data2.x1, data2.x2,seriestype=:scatter)
```

Now we write a function that will take a set of weights, update them and return the updated weights.
```
# A function to train weights. This is not the normal perceptron function but changes weights by an amount proportional to
# the distance fromthe seperation line. This ensures that the weights move more, the greater the seperation line is oway from
# correctly seperating a datapoint.
function IteratePerceptronWeights(w, LearningRate = 1)
    for i in 1:length(data[:y])
        target = data[i,:y]
        score = w[1] + (w[2]*data[i,:x1]) + (w[3]*data[i,:x2])
        ypred = 2*((score > 0)-0.5)
        if abs(target-ypred) > 1e-10
            update = LearningRate * -sign(score) * sqrt(abs(score))
            w[1] = w[1] + update
            w[2] = w[2] + update*data[i,:x1]
            w[3] = w[3] + update*data[i,:x2]
        end
    end
    return(w)
end
InitialGuess = [1.0, -2.0, 0.5]
FP = fixed_point(IteratePerceptronWeights, InitialGuess; Algorithm = MPE, PrintReports = true)
```

We can verify that the set of weights represented by the fixed\_point function does correctly seperate the data by plotting it:
```
# Plotting new seperation line
x1 = -6.0:0.1:6.0
w = FP.FixedPoint_
x2_on_sep_line = (-w[1] .- w[2] .* x1) ./ w[3]
plot!(x1,x2_on_sep_line, label ="SeperationLine")
```

## 4.3 A consumption smoothing problem

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

This takes 22 iterates with the anderson algorithm which is drastically better than the 459 iterates it takes with the simple method.
