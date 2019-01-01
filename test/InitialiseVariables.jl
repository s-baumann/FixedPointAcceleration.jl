using FixedPointAcceleration
#func(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
function func(x::Array{Float64,1})
    # The first coordinate convergences to 4.0 by 1 unit per iterate.
    output = Array{Float64,1}(x)
    if abs(x[1] - 4.0) <= 1.0
        output[1] = 4.0
    elseif x[1] > 4.0
        output[1] = x[1] - 1.0
    else
        output[1] = x[1] + 1.0
    end
    # The second does aitken convergence to 2.3
    output[2] = x[2] + (2.3-x[2])/2.0
    return output
end
Inputs = Array{Float64,2}(undef,2,1)
Inputs[:,1] = [19.0,10.0]
Outputs = Array{Float64,2}(undef,size(Inputs)[1],0)
Algorithm = SEA
ConvergenceMetric = supnorm(Resids::Array{Float64, 1}) = maximum(abs.(Resids))
ConvergenceMetricThreshold = 0.000000001
MaxIter = 1000
MaxM = 10
ExtrapolationPeriod = 7
Dampening = 1.0
PrintReports = false
ReportingSigFig = 5
ConditionNumberThreshold = 1e3
ReplaceInvalids = ReplaceElements
