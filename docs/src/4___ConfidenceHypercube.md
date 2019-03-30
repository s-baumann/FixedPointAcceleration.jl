# 4 Applications
## 4.5 Finding a confidence hypercube for a multivariate normal distribution.

We can find a confidence interval that includes x\% of a univariate normal distribution easily. It is more difficult however to come up with a confidence interval (or confidence area) for a multivariate Gaussian distribution. The first reason is that there is some ambiguity in what such a confidence area should look like. Considering some dimensions of the distribution will be correlated, it may be natural to look for an eliptically shaped area that will give the smallest possible area that covers x\% of the probability mass of a multivariate normal distribution (See Korpela et al. 2017).

Parameterising an elliptical area may be difficult however and it may be more natural to define a hypercube based on some basis of the multivariate normal distribution. A natural algorithm to do this would be to guess cutoff points marking the edges of the hypercube. Integrate the pdf of the normal distribution over this hypercube. If the integral deviates from that desired then come up with a new guess with different cutoff points. Iterate until there is a hypercube containing the desired x\% of the mass of the distribution.

```
# Generating data
using Distributions
using FixedPointAcceleration
using HCubature
using Random

# Generating an example distribution.
# Without loss of generality we use means of zero in every dimension.
# We use a random matrix sampled from the Wishart distribution.
Random.seed!(1234)
dist = Normal()
Endowments = rand(LogNormal(), G, N)
Tastes      = rand(G, N)

```
