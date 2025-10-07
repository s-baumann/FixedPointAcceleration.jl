"""
This function performs Minimal Polynomial extrapolation (MPE) or Reduced Rank Extrapolation (RRE) given a matrix of previous iterates of the function.
### Inputs
 * `Iterates` - A matrix of inputs and outputs excluding jumps. Can be pieced together from Inputs and Outputs matrices of the fixed_point function using the put_together_without_jumps function
 * `Algorithm` - The Algorithm for polynomial extrapolation. Should be either "MPE" for minimal polynomial extrapolation or "RRE" for reduced rank extrapolation.
### Returns
 * A `Vector` containing the extrapolated vector.
"""
function PolynomialExtrapolation(Iterates::AbstractArray{R,2}, Algorithm::Symbol) where R<:Number
    if (Algorithm == :MPE)
        TotalColumnsOfIterates = size(Iterates)[2]
        OldDifferences         = Iterates[:,2:(TotalColumnsOfIterates-1)] .- Iterates[:,1:(TotalColumnsOfIterates-2)]
        LastDifference         = Iterates[:,TotalColumnsOfIterates]       .- Iterates[:,(TotalColumnsOfIterates-1)]
        InverseOldDifferences  = pinv(OldDifferences)
        cVector                = -InverseOldDifferences * LastDifference
        cVector                = vcat(cVector,1)
        sumVec                 = sum(cVector)
        return (Iterates[:,2:TotalColumnsOfIterates] * cVector) ./ sumVec
    elseif (Algorithm == :RRE)
        TotalColumnsOfIterates   = size(Iterates)[2]
        FirstColumn              = Iterates[:,1]
        Differences              = Iterates[:,2:(TotalColumnsOfIterates)]      - Iterates[:,1:(TotalColumnsOfIterates-1)]
        SecondDifferences        = Differences[:,2:(TotalColumnsOfIterates-1)] - Differences[:,1:(TotalColumnsOfIterates-2)]
        FirstDifference          = Differences[:,1]
        Differences              = Differences[:,1:(TotalColumnsOfIterates-2)]
        InverseSecondDifferences = pinv(SecondDifferences)
        return FirstColumn - ((Differences * InverseSecondDifferences) * FirstDifference)
    else
        error("Invalid Algorithm input. PolynomialExtrapolation function can only take Algorithm as :MPE or :RRE.")
    end
end

"""
This is a helper function for EpsilonExtrapolation
### Inputs
 * `Iterates` - A matrix representing different iterates with one iterate per column. Can be pieced together from Inputs and Outputs matrices of the fixed_point function using the put_together_without_jumps function
 * `Algorithm` - Algorithm for epsilon extrapolation. Should be either "VEA" for the vector extrapolation algorithm or "SEA" for the scalar epsilon algorithm.
### Returns
 *  A `Vector` with the extrapolated vector.
"""
function EpsilonExtrapolation(Iterates::AbstractArray{R,2}, Algorithm::Symbol) where R<:Number
    # The function cannot do anything to a one column input so will return input unchanged.
    if (size(Iterates)[2] == 1) return Iterates end
    if (size(Iterates)[2] % 2 == 0) Iterates = Iterates[:,2:size(Iterates)[2]] end
    if (!(Algorithm in [:VEA, :SEA])) error("Invalid Algorithm input. EpsilonExtrapolation function can only take Algorithm as VEA or SEA") end
    Mat = Iterates
    RowsOfMatrix    = size(Mat)[1]
    TotalColumnsOfMatrix = size(Mat)[2]
    PreviousMatrix = zeros(RowsOfMatrix,(TotalColumnsOfMatrix-1))
    for MatrixColumn in reverse(2:TotalColumnsOfMatrix)
        DiffMatrix = Mat[:,2:MatrixColumn] .- Mat[:,1:(MatrixColumn-1)]
        NewMatrix = PreviousMatrix + EpsilonExtrapolationVectorOfInverses(DiffMatrix, Algorithm)
        PreviousMatrix = Mat[:,2:(MatrixColumn-1)]
        Mat = NewMatrix
    end
    # The function can get NAs from the inversion (ie if differenceMatrix contains a zero
    # for SEA then there is division by zero). To avert this we try with 2 less columns.
    if eltype(Mat) <: Complex
        if (any(isnan.(real.(Mat))) | any(isnan.(imag.(Mat))) | any(ismissing.(Mat)))
            Mat = EpsilonExtrapolation(Iterates[:,3:size(Iterates)[2]],Algorithm)
        end
    else
        if (any(isnan.(Mat)) | any(ismissing.(Mat)))
            Mat = EpsilonExtrapolation(Iterates[:,3:size(Iterates)[2]],Algorithm)
        end
    end
    return Mat[:,1]
end

"""
This is a helper function for EpsilonExtrapolation
### Inputs
 * `DifferenceMatrix` - The matrix of the differences in elements to be inverted.
 * `Algorithm` - SEA or VEA.
### Returns
 * A `Vector` of the result of inverting each (column) vector in a mmatrix.
"""
function EpsilonExtrapolationVectorOfInverses(DifferenceMatrix::AbstractArray{T,2}, Algorithm::Symbol) where T<:Number
    if (size(DifferenceMatrix)[1] == 1) | (Algorithm == :SEA)
        return 1 ./ DifferenceMatrix
    else
        invs = transpose(pinv(DifferenceMatrix[:,1]))
        if size(DifferenceMatrix)[2] < 2 return invs end
        for i in 2:size(DifferenceMatrix)[2]
            invs = hcat(invs, transpose(pinv(DifferenceMatrix[:,i])))
        end
        return invs
    end
end
