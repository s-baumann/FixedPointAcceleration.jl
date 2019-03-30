# 4 Applications
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
