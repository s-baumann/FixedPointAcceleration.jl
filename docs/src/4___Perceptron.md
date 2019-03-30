# 4 Applications
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
