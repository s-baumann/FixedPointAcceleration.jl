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

Most acceleration algorithms perform poorly in accelerating the convergence of this perceptron training algorithm. This is due to the perceptron often converging by a ﬁxed increment. This occurs because multiple iterates can result in the same observations being misclassiﬁed and hence the same changeintheweights. Asaresultwewillusethesimplemethodwhichisguaranteedtobeconvergent for this problem (Novikoff, 1963).

First we generate a dataset:
```
# Generating linearly separable data
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
function IteratePerceptronWeights(w, LearningRate = 1)
    for i in 1:length(data[:y])
        target = data[i,:y]
        score = w[1] + (w[2]*data[i,:x1]) + (w[3]*data[i,:x2])
        ypred = 2*((score > 0)-0.5)
        if abs(target-ypred) > 1e-10
            update = LearningRate * 0.5*(target-ypred)
            w[1] = w[1] + update
            w[2] = w[2] + update*data[i,:x1]
            w[3] = w[3] + update*data[i,:x2]
        end
    end
    return(w)
end
InitialGuess = [1.0, -2.0, 0.5]
FP = fixed_point(IteratePerceptronWeights, InitialGuess; Algorithm = Simple, PrintReports = true)
```

Only the simple method is convergent here and it is relatively slow taking 1121 iterations. We can still get a beneﬁt from accelerators however if we can modify the training algorithm to give training increments that change depending on distance from the ﬁxed point. This can be done by updating the weights by an amount proportional to a concave function of the norm of $$w_0 + \sum_{i=1}^N w_i x_{i,j}$$. Note that the instances in which the weights are not updated stay the same and hence the modiﬁed training function will result in the same set of ﬁxed points as the basic function. This is done in the next piece of code where the MPE method is used. It can be seen that there is a substantial increase in speed with only 54 iterations required by the MPE method.

```
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

[^7]: Note that for perceptrons there are always uncountably many such fixed points where the perceptron correctly classifies the entire training set and will not further update. On the other hand it is possible that the data is not linearly separable in which case there may be no fixed point and the weights will continue to update forever.

## 4.3 Expectation Maximisation

Consider we have a set of data which has come from two different multivariate normal distributions. There is a probability $$\tau$$ that a datapoint comes from the first multivariate distribution.
```
# Generating data from two two-dimensional gaussian processes
using Distributions
using FixedPointAcceleration
using Random
using DataFrames
true_tau = 0.6
nobs_1 = 400
nobs_2 = convert(Int, round(nobs_1 * ((1-true_tau)/true_tau)))
Random.seed!(1234)
mu_1 = [0.0,8.0]
cov_1 = [2.0,0.5,2.0]
covar_1 = Symmetric([cov_1[1] cov_1[2]; cov_1[2] cov_1[3]])
md_1 = MultivariateNormal(mu_1,covar_1)
mu_2 = [-4.0,10.0]
cov_2 = [2.0,-0.75,12.0]
covar_2 = Symmetric([cov_2[1] cov_2[2]; cov_2[2] cov_2[3]])
md_2 = MultivariateNormal(mu_2,covar_2)

rands_from_1 = transpose(rand(md_1, nobs_1))
rands_from_2 = transpose(rand(md_2, nobs_2))
data1 = DataFrame([rands_from_1[:,1], rands_from_1[:,2]], [:x1, :x2])
data2 = DataFrame([rands_from_2[:,1], rands_from_2[:,2]], [:x1, :x2])
dd  = vcat(data1,data2)
# Plotting it:
plot(data1.x1, data1.x2,seriestype=:scatter)
plot!(data2.x1, data2.x2,seriestype=:scatter)
```
Now we want to estimate the parameter $$\tau$$, the means (represented above by mu\_1 and mu\_2) and the covariance matrices (represented above by cov\_1, cov\_2) using only the realised datapoints in the DataFrame called dd. We will refer to these parameters as $$\theta$$.

If we knew from which distribution each datapoint came, the above task would be considerably easier. We could separate the data and for each use [standard techniques](https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices) to find the expected mean and covariance matrix. We do not know from which distribution each datapoint came from however. We could use a guess for $$\theta$$ to estimate the probabilities of each datapoint coming from each distribution however (and call this vector of estimates by $$Z$$). Then we could choose maximum likelihood estimates of $$\theta$$ using our estimates of $$Z$$ in the likelihood expression. This is the EM approach. Note that it lends itself well to fixed point acceleration - We can write a function that given $$\theta$$ creates estimated probabilities of source distribution for each datapoint ($$Z$$) and uses these in a maximum likelihood expression to improve the estimate of $$\theta$$.

In this multivariate gaussian case there are simple expressions to choose parameters to maximise the likelihood. These are recounted on the [wikipedia article on expectation maximisation](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm) and are used below:
```
function z_estimate_given_theta(x::Array{Float64,1}, md_1::MultivariateNormal, md_2::MultivariateNormal, tau::Float64)
    pdf_1 = pdf(md_1, x)
    pdf_2 = pdf(md_2, x)
    return tau*pdf_1 / (tau*pdf_1 + (1-tau)*pdf_2)
end

function update_tau(Z::Array{Float64,1})
    return mean(Z)
end

function update_mu(dd::DataFrame, Z::Array{Float64,1})
    X = convert(Array{Float64,2}, dd[[:x1, :x2]])
    sum_Z = sum(Z)
    updated_mu = (transpose(Z) * X) ./sum_Z
    return vec(updated_mu)
end

function update_cov(dd::DataFrame, updated_mu::Array{Float64,1}, Z::Array{Float64,1})
    X_minus_mu = convert(Array{Float64,2}, dd[[:x1, :x2]]) .- transpose(updated_mu)
    sum_Z = sum(Z)
    updated_cov = (transpose(Z .* X_minus_mu) * X_minus_mu) ./sum_Z
    return [updated_cov[1,1], updated_cov[1,2], updated_cov[2,2]]
end

function update_theta(theta::Array{Float64,1}, dd::DataFrame)
    # We will use the convention that theta's 11 entries are (mu_1, cov_1, mu_2, cov_2, tau). First unpacking theta:
    mu_1    = theta[[1,2]]
    cov_1   = theta[[3,4,5]]
    covar_1 = Symmetric([cov_1[1] cov_1[2]; cov_1[2] cov_1[3]])
    md_1 = MultivariateNormal(mu_1,covar_1)
    mu_2    = theta[[6,7]]
    cov_2   = theta[[8,9,10]]
    covar_2 = Symmetric([cov_2[1] cov_2[2]; cov_2[2] cov_2[3]])
    md_2 = MultivariateNormal(mu_2,covar_2)
    tau     = theta[11]
    # Getting Z
    Z = Array{Float64,1}(undef,size(dd)[1])
    for i in 1:size(dd)[1]
        Z[i] = z_estimate_given_theta([dd[i,:x1], dd[i,:x2]], md_1, md_2, tau)
    end

    # Updating Tau
    updated_tau = update_tau(Z)
    # Updating mu1
    updated_mu_1 = update_mu(dd,Z)
    updated_mu_2 = update_mu(dd, 1 .- Z)
    # Updating Cov
    updated_cov_1 = update_cov(dd, updated_mu_1, Z)
    updated_cov_2 = update_cov(dd, updated_mu_2, 1 .- Z)
    # Returning theta
    updated_theta = vcat(updated_mu_1, updated_cov_1, updated_mu_2, updated_cov_2, updated_tau)
    return updated_theta
end
```
Now we can come up with a choice for an initial guess based on eyeballing the plotted data. We can then put it into the fixed\_point function to get ML estimates of these distributional parameters as well as $$\tau$$:

```
InitialGuess = [0.5, 7.5, 2.0, 0.0, 2.0, -5.0, 7.5, 2.0, 0.0, 10.0, 0.5]
fp_anderson = fixed_point(x -> update_theta(x,dd), InitialGuess; Algorithm = Anderson, PrintReports = true)
fp_simple   = fixed_point(x -> update_theta(x,dd), InitialGuess; Algorithm = Simple, PrintReports = true)
```
We can see that the Anderson method only takes 15 iterations while the simple method takes 80. By checking the generated fixedpoint against the data generation process it can also be verified that the fixedpoint we find provides quite good estimates.

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

## 4.5 Finding a confidence hypercube for a multivariate normal distribution.

We can find a confidence interval that includes x% of a univariate normal distribution easily. It is more difficult however to come up with a confidence interval (or confidence area) for a multivariate Gaussian distribution. The first reason is that there is some ambiguity in what such a confidence area should look like. Considering some dimensions of the distribution will be correlated, it may be natural to look for an elliptically shaped area that will give the smallest possible area that covers x% of the probability mass of a multivariate normal distribution (See Korpela et al. 2017).

Parameterising an elliptical area may be difficult however and it may be more natural to define a hypercube based on some basis of the multivariate normal distribution. A natural algorithm to do this would be to:
1. Guess cutoff points marking the edges of the hypercube.
2. Integrate the pdf of the normal distribution over this hypercube.
3. If the integral deviates from that desired then come up with a new guess with different cutoff points.
Then we iterate these steps until there is a hypercube containing the desired x\% of the mass of the distribution. While this should reach a fixed point of cutoff points there will be many such hypercubes. For instance we could take the x% confidence interval off the marginal distribution (in one dimension) of the multivariate normal. Using these cutoffs for this dimension and $[-\infty, \infty]$ as cutoffs for all of the other dimensions we will have a x% confidence area.

We are likely to be more interested in the most "central" hypercube. We shall put in additional restriction that in each dimension the hypercube should extend the same number of standard deviations above and below the mean.

First we generate an example multivariate normal distribution:
```
using Distributions
using FixedPointAcceleration
using HCubature
using Random
using LinearAlgebra

# Generating an example distribution.
# Without loss of generality we use means of zero in every dimension.
# We use a random matrix sampled from the Wishart distribution.
Random.seed!(1234)
prob_means = repeat([0.0],100)
dims = 100
wish = Wishart(dims, diagm(0 => ones(dims)))
covar_matrix = rand(wish, 1)[1]
dist = MvNormal(prob_means, Symmetric(covar_matrix))
chol_of_covar_matrix = LowerTriangular(cholesky(covar_matrix).L)
```

Now normal numerical integration routines does not scale well with the number of dimensions and this distributions has 100 dimensions. As a result we will instead infer the integral through the use of the Sobol sequence:
```
# We create an array of values for integrating cheaply
using Sobol
function get_sobol_draws(chol, num::Int, sob_seq::SobolSeq)
    dims = size(chol)[1]
    array = Array{Float64,2}(undef, num, dims)
    for i in 1:num
        sobs = next!(sob_seq)
        normal_draw = quantile.(Ref(Normal()), sobs)
        scaled_draw = chol * normal_draw
        array[i,:] = scaled_draw
    end
    return array
end
draws = get_sobol_draws(chol, 100000, SobolSeq(dims))
```

We can now write a function that updates the number of standard deviations above and below the mean that defines the edges of the hypercube.
```
# Update function
function one_iterate(cutoff_multiplier::Float64, target::Float64; tuning_parameter::Float64 = 1.0)
    cutoffs = cutoff_multiplier .* sqrt.(diag(covar_matrix))
    number_of_draws = size(draws)[1]
    in_confidence_area = 0
    for i in 1:number_of_draws
        in_confidence_area += all(abs.(draws[i,:]) .< cutoffs)
    end
    mass_in_area = in_confidence_area/number_of_draws
    confidence_gap = target - mass_in_area
    return cutoff_multiplier + confidence_gap * tuning_parameter
end
FP = fixed_point(x -> one_iterate.(x, 0.95), [2.0]; Algorithm = Anderson, PrintReports = true)
# The final number of standard deviations above/below the mean to use is stored in FP:
cutoff_multiplier = FP.FixedPoint_[1]
# We can find the upper and lower edges of the hypercube in each dimension. They are stored in each dimension in the below array of tuples.
cutoffs = vcat(zip(-cutoff_multiplier .* sqrt.(diag(covar_matrix)) , cutoff_multiplier .* sqrt.(diag(covar_matrix)))...)
```

## 4.6 Importance Sampling

To lower the variance of a Monte Carlo estimation, importance sampling is often used. We can use a simple method of importance sampling for function integrated over a multivariate standard normal distribution using the algorithm described in Glasserman (2003) on page 268 of that book.

To briefly recount it we want to find $E[G(x)]$ where $x \in \Re^d$ is $N(0,I)$. We have the change of measure to a measure called $\mu$ which:

$$E[G(Z)] &= E_\mu \left[ G(Z)e^{-\mu^\prime Z + \frac{1}{2}\mu^\prime\mu} \right]$$

for any $\mu \in \Re^d$ where $d$ is the dimensionality of the function $G(\cdot)$. We can simulate this with the algorithm:

for paths $i = 1, ..., N$\\
	 generate $Z_i \sim N(\mu, I)$\\
	 $Y_i \leftarrow G(Z_i) \exp(-\mu^\prime Z_i + \frac{1}{2}\mu^\prime \mu)$\\
return $\frac{\sum_{i=1}^N Y_i}{N}$

Now we need to figure out the vector $\mu$ which is composed of the shifts in the mean for each normal variable and thereby represents the change in probability measure. Note that for any vector our estimator should be unbiased and consistent but some can be more efficient than others. In the special case\footnote{Which is satisfied for CVA and DVA but not for FVA, thus the scope of this paper.} where $G(x) \geq 0 \forall x \in \Re^d$, an efficient choices is the vector $\mu$ which makes the following equation hold (for working out see equation 4.89 of Glasserman):

$$\Delta F(\mu) = \mu$$

Where $F(x) = \ln(G(x))$ and  $\Delta F(\mu)$ is the Jacobian of the function $F(\cdot)$

The first problem here is how to to efficiently get the Jacobian. For a high dimensional problem (which all XVA problems become) numerical differentiation will not work but aad will work so we can use the ForwardDiff package. The second problem is the high dimensional fixedpoint problem for which we can used FixedPointAcceleration.

```
using LinearAlgebra
using Distributions
using FixedPointAcceleration
using ForwardDiff
using Random
twister = MersenneTwister(1)
dims = 5
random_function_multiples = rand(twister, dims)
Wish = Wishart(dims, diagm(0 => ones(dims)))
random_PSD_matrix = rand(twister, Wish)
chol = cholesky(random_PSD_matrix).factors

function G(x::Array)
    transformed_normals = 0.3 .* (chol * x)
    return sum(random_function_multiples .* exp.(transformed_normals))
end
function F(x::Array)
    return log(G(x))
end
function grad_F(x::Array)
  return ForwardDiff.gradient(F, x)
end

fp = fixed_point(grad_F, repeat([1.0],dims))
drifts = fp.FixedPoint_
iterations = fp.Iterations_

batches = 200000
batch_size = 30
paths = batches * batch_size
draws = rand(twister, Normal(), paths, dims)

normal_batch_result = Array{Float64,1}(undef, batches)
is__batch_result     = Array{Float64,1}(undef, batches)
normal_result = Array{Float64,1}(undef, batches)
is_result     = Array{Float64,1}(undef, batches)

for batch in 1:batches
    println("Now doing batch ", batch)
    batch_vec = Array{Float64,1}(undef, batch_size)
    batch_is_vec = Array{Float64,1}(undef, batch_size)
    for i in 1:batch_size
        row = (batch-1)*batch_size + i
        vv = draws[row,:]
        batch_vec[i] = G(vv)
        vv2 = vv + drifts
        batch_is_vec[i] = G(vv2) * exp(0.5 * (drifts' * drifts) - drifts' * vv2)
    end
    normal_batch_result[batch] = mean(batch_vec)
    is__batch_result[batch]    = mean(batch_is_vec)
    normal_result[batch]       = mean(normal_batch_result[1:batch])
    is_result[batch]           = mean(is__batch_result[1:batch])
end

final_val = normal_result[batches]
upto = 3000

using Plots
plt = plot(1:upto , normal_result[1:upto] , label= "Normal QMC")
plt = plot!(iterations .+ collect(1:upto), is_result[1:upto], label= "QMC with Importance Sampling")
plt = plot!(1:(upto+iterations), repeat([final_val], upto+iterations), label= "Converged Value")
```

You can see that convergence is much faster here using importance sampling.
