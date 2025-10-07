"""
This function takes the previous inputs and outputs and assembles a matrix with both excluding jumps.
### Inputs
 * `Inputs` - This is an N x A matrix of previous inputs for which corresponding outputs are available. In this case N is the dimensionality of the fixed point vector that is being sought (and each column is a matrix that is input to the "Function") and A is the number of previous Inputs/Outputs that are being provided to the fixed point.
 * `Outputs` This is a matrix of "Function" values for each column of the "Inputs" matrix.
 * `AgreementThreshold` - A parameter for determining when a column in Inputs and a column in Outputs match. They are deemed to match if the sum of the absolute values of the difference in the columns is less than AgreementThreshold.
### Returns
 * A matrix of inputs and outputs excluding jumps.
"""
function put_together_without_jumps(Inputs::AbstractArray{T,2}, Outputs::AbstractArray{T,2}, AgreementThreshold::Float64 = 1e-10) where T<:Real
  if (any(size(Inputs) != size(Outputs))) error("Inputs and Outputs matrices are not comformable.") end
  size_of_dims = size(Inputs)

  if (size_of_dims[2] == 1) return(hcat(Inputs, Outputs)) end
  Difference = (Inputs[:,2:(size_of_dims[2])] .- Outputs[:,1:(size_of_dims[2]-1)])
  Sum_Of_Differences = sum(Difference, dims = 1)[1,:]
  Agreements = abs.(Sum_Of_Differences) .< AgreementThreshold
  if (all(Agreements))
      return hcat(Inputs[:,1:(size_of_dims[2])], Outputs[:, size_of_dims[2]])
  else
      LocationsOfBreaks = findall(Agreements .== false)
      LastBreak = LocationsOfBreaks[length(LocationsOfBreaks)]
      return hcat(Inputs[:,(LastBreak+1):(size_of_dims[2])] , Outputs[:, size_of_dims[2]])
  end
end
