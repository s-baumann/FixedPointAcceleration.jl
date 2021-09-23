using FixedPointAcceleration
using Test

func(x) = cos.(x)
Inputs = 1.1
fp_simple   = fixed_point(func, Inputs; Algorithm = :Simple)
fp_simple.Convergence_ < 1e-10
fp_anderson = fixed_point(func, Inputs; Algorithm = :Anderson)
fp_anderson.Convergence_ < 1e-10
fp_aitken   = fixed_point(func, Inputs; Algorithm = :Aitken)
fp_aitken.Convergence_ < 1e-10
fp_newton   = fixed_point(func, Inputs; Algorithm = :Newton)
fp_newton.Convergence_ < 1e-10
fp_VEA      = fixed_point(func, Inputs; Algorithm = :VEA)
fp_VEA.Convergence_ < 1e-10
fp_SEA      = fixed_point(func, Inputs; Algorithm = :SEA)
fp_SEA.Convergence_ < 1e-10
fp_MPE      = fixed_point(func, Inputs; Algorithm = :MPE)
fp_MPE.Convergence_ < 1e-10
fp_RRE      = fixed_point(func, Inputs; Algorithm = :RRE)
fp_RRE.Convergence_ < 1e-10

# Testing input of inputs without outputs.
@test_logs (:warn,"If you do not give outputs to the function then you can only give one vector of inputs (in a 2d array) to the fixed_pointFunction. So for a function that takes an N dimensional array you should input a Array{Float64}(N,1) array.  As you have input an array of size Array{Float64}(N,k) with k > 1 we have discarded everything but the last column to turn it into a Array{Float64}(N,1) array.\n") fixed_point(func, fp_RRE.Inputs_; Algorithm = :RRE)

# Testing input of Outputs with different shape to inputs.
@test_logs (:warn,"If you input a matrix of outputs as well as a matrix of inputs then inputs and outputs must be the same shape. As they differ in this case the last column of the inputs matrix has been taken as the starting point and everything else discarded.") fixed_point(func, fp_RRE.Inputs_; Outputs = fp_RRE.Outputs_[:,2:3], Algorithm = :RRE)
