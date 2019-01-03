using SchumakerSpline
using HCubature
using Distributions
using Random
using Optim
using FixedPointAcceleration
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

function new_point(level::Float64, stepp::Float64, max_grad::Float64, grad_frac::Float64)
    expp = exp(grad_frac)
    frac = expp / (expp+1)
    if isinf(max_grad)
        return level + expp
    else
        return level + stepp * max_grad * frac
    end
end

function invert_new_point(level::Float64, previous_level::Float64, previous_stepp::Float64, previous_max_grad::Float64)
    rise = level - previous_level
    if isinf(previous_max_grad)
        return log(rise)
    else
        R = rise/(previous_stepp * previous_max_grad)
        return log(R/(1-R))
    end
end

function ShapeToBud(ShapeVec::Array{Float64,1})
    lenlen = length(ShapeVec)
    BudgetValues = Array{Float64,1}(undef, lenlen)
    if lenlen > 0 BudgetValues[1] = ShapeVec[1] end
    if lenlen > 1 BudgetValues[2] = new_point(BudgetValues[1], 1.0,Inf,ShapeVec[2]) end
    if lenlen > 2
        for i in 3:lenlen
            stepp = BudgetStateSpace[i] - BudgetStateSpace[i-1]
            previous_grad = (BudgetValues[i-1]-BudgetValues[i-2])/(BudgetStateSpace[i-1] - BudgetStateSpace[i-2])
            BudgetValues[i] = new_point(BudgetValues[i-1], stepp, previous_grad, ShapeVec[i])
        end
    end
    return BudgetValues
end

function BudToShape(BudgetValues::Array{Float64,1})
    lenlen = length(BudgetValues)
    new_shape_vec = Array{Float64,1}(undef, lenlen)
    if lenlen > 0 new_shape_vec[1] = BudgetValues[1] end
    if lenlen > 1 new_shape_vec[2] = invert_new_point(BudgetValues[2],BudgetValues[1], 1.0,Inf) end
    if lenlen > 2
        for i in 3:lenlen
            stepp = BudgetStateSpace[i] - BudgetStateSpace[i-1]
            previous_grad = (BudgetValues[i-1]-BudgetValues[i-2])/(BudgetStateSpace[i-1] - BudgetStateSpace[i-2])
            new_shape_vec[i] = invert_new_point(BudgetValues[i],BudgetValues[i-1], stepp, previous_grad)
        end
    end
    return new_shape_vec
end

function Reparameterised_FP_Vector(ShapeVec::Array{Float64,1})
    BudgetValues = ShapeToBud(ShapeVec)
    new_Budget_Values = OneIterateBudgetValues(BudgetValues)
    new_shape_vec = BudToShape(new_Budget_Values)
    return new_shape_vec
end
function shapeconvergence(Inputs, Outputs)
    return maximum(abs.(ShapeToBud(Inputs) .- ShapeToBud(Outputs)))
end

# Testing that the conversions work properly.
fp = fixed_point(OneIterateBudgetValues, InitialGuess; PrintReports = true, MaxIter = 1)
shape_guess = BudToShape(InitialGuess)
fp_reparam = fixed_point(Reparameterised_FP_Vector, shape_guess; PrintReports = true, MaxIter = 1, ConvergenceMetric = shapeconvergence)
iterated = ShapeToBud(fp_reparam.Outputs_[:,1])
sum(abs.(fp.Outputs_[:,1] .- iterated) .> 1e-10)  == 0

# fixed point acceleration with the reparameterised version.
fp = fixed_point(OneIterateBudgetValues, InitialGuess; PrintReports = true, ConvergenceMetricThreshold = 1e-06)
fpfp = fp.FixedPoint_
fp_ander = fixed_point(Reparameterised_FP_Vector, shape_guess; PrintReports = true, ReportingSigFig = 10, ConvergenceMetricThreshold = 1e-06, ConvergenceMetric = shapeconvergence)
iterated = ShapeToBud(fp_ander.FixedPoint_)
sum(abs.(fpfp .- iterated) .> 1e-4)  == 0
fp_aitken = fixed_point(Reparameterised_FP_Vector, shape_guess; PrintReports = true, ReportingSigFig = 10, ConvergenceMetricThreshold = 1e-06, Algorithm = Aitken, ConvergenceMetric = shapeconvergence)
iterated = ShapeToBud(fp_aitken.FixedPoint_)
sum(abs.(fpfp .- iterated) .> 1e-4)  == 0
fp_newton = fixed_point(Reparameterised_FP_Vector, shape_guess; PrintReports = true, ReportingSigFig = 10, ConvergenceMetricThreshold = 1e-06, Algorithm = FixedPointAcceleration.Newton, ConvergenceMetric = shapeconvergence)
iterated = ShapeToBud(fp_newton.FixedPoint_)
sum(abs.(fpfp .- iterated) .> 1e-4)  == 0
fp_sea = fixed_point(Reparameterised_FP_Vector, shape_guess; PrintReports = true, ReportingSigFig = 10, ConvergenceMetricThreshold = 1e-06, Algorithm = SEA,  ConvergenceMetric = shapeconvergence)
iterated = ShapeToBud(fp_sea.FixedPoint_)
sum(abs.(fpfp .- iterated) .> 1e-4)  == 0
