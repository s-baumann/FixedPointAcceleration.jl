

func(x) = sqrt(x)
Inputs = ones(7,1)*16
Inputs = hcat(Inputs, ones(7,1)*4)
Outputs = hcat(ones(7,1)*4, ones(7,1)*2)
Algorithm = :Simple
ConvergenceMetric(Resids::Array{Float64, 1}) = maximum(abs.(Resids))
ConvergenceMetricThreshold = 1e-10
MaxIter = 1e3
MaxM = 10
ExtrapolationPeriod = 7
Dampening = 1.0
PrintReports = false
ReportingSigFig = 5
ConditionNumberThreshold = 1e3
Plot = :NoPlot
ConvergenceFigLags = 5
ChangePerIteratexaxis = []

#function aa( Inputs::Array{Float64, 2}, Outputs::Array{Float64, 2} )
#Inputs[:,1] + Outputs
#end
#aa(Inputs, Outputs)
x = [-1.0,0.0,1.0]
