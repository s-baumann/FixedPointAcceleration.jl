using FixedPointAcceleration

Function(x) = sqrt(x)
Inputs = ones(7,1)*16
Inputs = hcat(Inputs, ones(7,1)*4)
Outputs = hcat(ones(7,1)*4, ones(7,1)*2)
Algorithm = ("Anderson")
ConvergenceMetric(Resids::Array{Float64}) = maximum(abs.(Resids))
ConvergenceMetricThreshold = 1e-10
MaxIter = 1e3
MaxM = 10
ExtrapolationPeriod = 7
Dampening = 1.0
PrintReports = false
ReportingSigFig = 5
ConditionNumberThreshold = 1e3
Plot = ("NoPlot", "ConvergenceFig", "ChangePerIterate")
ConvergenceFigLags = 5
ChangePerIteratexaxis = []
