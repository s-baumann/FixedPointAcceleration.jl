using FixedPointAcceleration
#func(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
func = simple_vector_function(x) = [0.5*sqrt(x[1] + x[2]), 1.5*x[1] + 0.5*x[2]]
Inputs = Array{Float64,2}(undef,2,1)
Inputs[:,1] = [0.3,900]
Outputs = Array{Float64,2}(undef,size(Inputs)[1],0)
Algorithm = Simple
ConvergenceMetric = supnorm(input::Array{Float64, 1}, output::Array{Float64,1}) = maximum(abs.(output .- input))
ConvergenceMetricThreshold = 0.000000001
MaxIter = 1000
MaxM = 10
ExtrapolationPeriod = 7
Dampening = 1.0
PrintReports = false
ReportingSigFig = 10
ConditionNumberThreshold = 1e3
ReplaceInvalids = ReplaceElements
