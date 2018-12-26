
func(x) = sqrt.(x)
Inputs = Array{Float64,1}([2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 25.6])
Outputs = func(Inputs)
Inputs = hcat(Inputs, Outputs[:,(2-1)])
result = func.(Inputs[:,2])
Outputs = hcat(Outputs, result)
Inputs = hcat(Inputs, Outputs[:,(3-1)])
result = func.(Inputs[:,3])
Outputs = hcat(Outputs, result)
Inputs = hcat(Inputs, Outputs[:,(4-1)])
result = func.(Inputs[:,4])
Outputs = hcat(Outputs, result)

Iterates = put_together_without_jumps(Inputs, Outputs)
Iterates == hcat(Inputs, Outputs[:,size(Outputs)[2]])

NewGuess = [5.4, 5.5, 2.3, 3.4, 1.2, 0.56, 9.0]
Inputs = hcat(Inputs, NewGuess)
result = func.(Inputs[:,5])
Outputs = hcat(Outputs, result)
Inputs = hcat(Inputs,  Outputs[:,(6-1)])
result = func.(Inputs[:,6])
Outputs = hcat(Outputs, result)

Iterates = put_together_without_jumps(Inputs, Outputs)
Iterates == hcat(NewGuess, func.(NewGuess), func.(func.(NewGuess)))
