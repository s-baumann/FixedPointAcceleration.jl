using FixedPointAcceleration
#func(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
func(x) = cos.(x)
Inputs = Array{Float64,2}(undef,1,1)
Inputs[1,1] = 1.1
Outputs = Array{Float64,2}(undef,size(Inputs)[1],0)
Algorithm = RRE
ConvergenceMetric = supnorm(Resids::Array{Float64, 1}) = maximum(abs.(Resids))
ConvergenceMetricThreshold = 0.000000001
MaxIter = 1000
MaxM = 10
ExtrapolationPeriod = 7
Dampening = 1.0
PrintReports = false
ReportingSigFig = 5
ConditionNumberThreshold = 1e3
ReplaceInvalids = NoAction
