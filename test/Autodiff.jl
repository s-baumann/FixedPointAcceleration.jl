using ForwardDiff
using FixedPointAcceleration
func(x, params) = [params[1]*sqrt(abs(x[1] + x[2])), 1.5*x[1] + params[2]*x[2]]

function ML_estimation(params)
    Inputs = [0.3,900.0]
    fp_simple   = fixed_point(inp -> func(inp, params), Inputs; Algorithm = Simple)
    return fp_simple.FixedPoint_[1]
end

params = [1.5, 0.5]
ML_estimation(params)
grads = ForwardDiff.gradient(ML_estimation, params)
isa(grads, Array{Float64,1})
