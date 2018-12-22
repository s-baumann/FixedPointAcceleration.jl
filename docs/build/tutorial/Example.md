
<a id='Example:-A-simple-consumption-smoothing-problem-1'></a>

## Example: A simple consumption smoothing problem


Consider a consumer that has a budget of $B_t$ at time $t$ and a periodic income of $1$. She has a periodic utility function given by:


$u_t = \epsilon_t x_t^{0.2}$

where $x_t$ is spending in period $t$ and $\epsilon_t$ is the shock in period $t$ drawn from some stationary nonnegative shock process with pdf $f(\epsilon)$.


The problem for the consumer in period $t$ is:


$V(B_t | \epsilon_{t}) =  \max_{0 < x_t < B} \hspace{0.5cm} \epsilon_t x_t^{0.2} + \beta E_t[ V(B_{t+1})]$

Where $\beta$ is a discounting factor and $B_{t+1} = 1 + B_t - x_t$.


<a id='Algorithm-1'></a>

### Algorithm


We can first note that due to the shock process it is not possible to get analytical equations to describe the paths of spending and the budget over the long term. We can get a numerical solution however. The key step is to find expressions for the expected value function as a function of $B_{t+1}$. With this we can run simulations with random shocks to see the long term distributions of $x_t$ and $B_t$. The algorithm we will use is:


1. We discretize the budget statespace.
2. We make an initial guess of the future value function $E[ V(B_{t+1})]$ at every value of $B_{t+1}$ in our discretized statespace. A deterministic approximation of the problem (assume shocks will be $E[\epsilon_{t}]$ forever) is often good for this.
3. We use the schumaker spline to join together our estimates at each point in the discretized statespace into an interpolation function.
4. At every point in the statespace we create updated estimates   * We use a cubature package to estimate $\int^\infty_{0} V(B_t | \epsilon_t ) f(\epsilon_t) d\epsilon_t$   * In this integral we use a one dimensional optimizer along with the above equation to evaluate $V(B_t | \epsilon_t )$
5. Check your convergence criteria. Are the new $V(B_t | \epsilon_t )$ values very close to the old values?   * If they are we end the algorithm   * If they are not we go back to point 3.  This strategy relies on the consumption problem being a contraction mapping. This means that if we use this algorithm we will converge to a fixed point function.  ### Why a schumaker spline is necessary  There are a few reasons we need the spline to be shape preserving and without optimization to make this work:    * Shape preservation is necessary so the spline is globally concave (the $V(B_t | \epsilon_t )$ values will always be concave as the periodic utility function is always concave). This is necessary for the one dimensional optimisation step. The periodic utility function is concave and the future value function needs to be concave (in $B_{t+1}$) to guarantee a unique local optima. Other splines can incorporate tiny convexities in the intervals between interpolation points which can lead to multiple local optima.      * There do exist shape preserving splines that incorporate an optimization step. These are not effective for this kind of problem because it is not possible to converge to a fixed point. In each iteration the optimization settles on a slightly different parameter set which means the future value function can "jump around". In addition the optimization step can take a significant amount of time because evaluations of the future value function spline take longer.

