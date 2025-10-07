"""
    fixed_point_new_input(Inputs::AbstractArray{T,2}, Outputs::AbstractArray{T,2}, Algorithm::Symbol = :Anderson;
                               MaxM::Integer = 10, SimpleStartIndex::Integer = 1, ExtrapolationPeriod::Integer = 1, Dampening::S = AbstractFloat(1), Dampening_With_Input::Bool = false,
                               ConditionNumberThreshold::R = AbstractFloat(1000), PrintReports::Bool = false, ReplaceInvalids::InvalidReplacement = :NoAction) where R<:Real where S<:Real where T<:Real

This function takes the previous inputs and outputs from the fixed_point function and determines what vector to try next in seeking a fixed point.
### Inputs
* `Inputs` - This is an N x A matrix of previous inputs for which corresponding outputs are available. In this case N is the dimensionality of the fixed point vector that is being sought (Hence each column is a matrix that is input to the "Function") and A is the number of previous Inputs/Outputs that are being provided to the fixed point.
* `Outputs` - This is a matrix of Function values for the each column of the `Inputs` matrix.
* `Algorithm` - This is the fixed point Algorithm to be used. It can be "Anderson", "Simple", "Aitken", "Newton", "MPE", "RRE", "VEA", "SEA".
* `MaxM` - This is the number of saved iterates that are used in the Anderson algorithm. This has no role if another Algorithm is used.
* `SimpleStartIndex` - This is the index for what column in the input/output matrices did the algorithm start doing simple iterates without jumps. This is used for all Algorithms except the simple and Anderson Algorithms where it has no effect.
* `ExtrapolationPeriod` - This is the number of simple iterates to perform before extrapolating. This is used for the MPE, RRE, VEA and SEA Algorithms and has no effect if another Algorithm is chosen.
* `Dampening` - This is the dampening parameter. By default it is 1 which means no dampening takes place. It can also be less than 1 (indicating dampening) or more than 1 (indicating extrapolation).
* `Dampening_With_Input` - This is a boolean that indicates whether the dampening parameter should be multiplied by the input (if true) or the output of the most recent iterate.
* `ConditionNumberThreshold` - This is what threshold should be chosen to drop previous iterates if the matrix is ill conditioned. Only used in Anderson acceleration.
* `PrintReports` - This is a boolean describing whether to print ongoing ConvergenceMetric values for each iterate.
### Returns
 * A `Vector` of the next guess for the fixed point.
### Examples
    FPFunction = function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])}
    A = fixed_point( Function = FPFunction, Inputs = [0.3,900], MaxIter = 6, Algorithm = :Simple)
    NewGuessAnderson = fixed_point_new_input(A[:Inputs], A[:Outputs], Algorithm = :Anderson)
    NewGuessVEA = fixed_point_new_input(A[:Inputs], A[:Outputs], Algorithm = :VEA)
    NewGuessMPE = fixed_point_new_input(A[:Inputs], A[:Outputs], Algorithm = :MPE)
    NewGuessAitken = fixed_point_new_input(A[:Inputs], A[:Outputs], Algorithm = :Aitken)
"""
function fixed_point_new_input(Inputs::AbstractArray{T,2}, Outputs::AbstractArray{T,2}, Algorithm::Symbol = :Anderson;
                               MaxM::Integer = 10, SimpleStartIndex::Integer = 1, ExtrapolationPeriod::Integer = 1, Dampening::S = 1.0, Dampening_With_Input::Bool = false,
                               ConditionNumberThreshold::R = 1000.0, PrintReports::Bool = false, ReplaceInvalids::Symbol = :NoAction) where R<:Real where S<:Number where T<:Number
    CompletedIters = size(Outputs)[2]
    simple_iterate = Outputs[:,CompletedIters]
    proposed_input = repeat([NaN], size(simple_iterate)[1])
    if Algorithm == :Simple
         proposed_input = simple_iterate
    elseif Algorithm == :Anderson
        if CompletedIters < 2
            if (PrintReports) print("                           Used:",  lpad(0, 3)," lags. ") end
            proposed_input = simple_iterate
        else
            VectorLength   = size(Outputs)[1]
            M = min(MaxM-1,CompletedIters-1,VectorLength)

            recent_Outputs  = Outputs[:, (CompletedIters-M):CompletedIters]
            recent_Inputs   = Inputs[ :, (CompletedIters-M):CompletedIters]
            Resid           = recent_Outputs .-  recent_Inputs
            DeltaOutputs    = recent_Outputs[:,2:(M+1)] .- recent_Outputs[:,1:M]
            DeltaResids     = Resid[:,2:(M+1)]   .- Resid[:,1:M]
            LastResid       = Resid[:,M+1]
            LastOutput      = recent_Outputs[:,M+1]
            Coeffs          = repeat([NaN], size(DeltaOutputs)[2])
            ConditionNumber = NaN
            while any(isnan.(Coeffs))
                if isempty(DeltaResids)
                    # This happens if there is convergence by constant increments and thus the
                    #  most recent DeltaResids is all zeros. So we end up dropping all DeltaResids and get an error here.
                    break
                end
                ConditionNumber = cond(DeltaResids)
                if ConditionNumber > ConditionNumberThreshold
                    M = M-1
                    DeltaOutputs= DeltaOutputs[:, 2:(M+1)]
                    DeltaResids = DeltaResids[ :, 2:(M+1)]
                    Coeffs      = repeat([NaN], size(DeltaOutputs)[2])
                    continue
                end
                # Handle complex numbers by using pinv instead of GLM.fit
                if eltype(DeltaResids) <: Complex
                    Coeffs = pinv(DeltaResids) * LastResid
                else
                    Fit = fit(LinearModel,  hcat(DeltaResids), LastResid)
                    Coeffs = Fit.pp.beta0
                end
                if any(isnan.(Coeffs))
                    M = M-1
                    if (M < 1.5)
                        # This happens occasionally in test cases where the iteration is very close to a fixed point.
                        if (PrintReports) print("                          Used:",  lpad(0, 3)," lags. ") end
                        break
                    end
                    DeltaOutputs = DeltaOutputs[:, 2:(M+1)]
                    DeltaResids  = DeltaResids[ :, 2:(M+1)]
                end
            end
            if isempty(Coeffs)
                if (PrintReports) print("Condition number is ", lpad("NaN", 5),". Used:",  lpad(0, 3)," lags. ") end
                proposed_input = repeat([NaN], VectorLength)
            else
                if (PrintReports) print("Condition number is ", lpad(round(ConditionNumber, sigdigits = 2), 5),". Used:",  lpad(M+1, 3)," lags. ") end
                proposed_input = LastOutput .- (Dampening .* vec(DeltaOutputs * Coeffs))
            end
        end
    elseif Algorithm == :Aitken
        if ((CompletedIters + SimpleStartIndex) % 3) == 0
            # If we are in 3rd, 6th, 9th, 12th iterate from when we started Acceleration then we want to do a jumped Iterate,
            # First we extract the guess that started this run of 3 iterates (x), the Function applied to it (fx) and the function applied to that (ffx)
            x = Inputs[: ,(CompletedIters -1)]
            fx = Outputs[:,(CompletedIters -1)]
            ffx = Outputs[:,CompletedIters]
            # Now using the appropriate formula to make a new guess. Note that if a vector is input here it is used elementwise.
            proposed_input = x .- ((fx .- x).^2 ./ (ffx .- 2 .* fx .+ x))
        else
            # We just do a simple iterate. We do an attempt with the latest iterate.
            proposed_input = simple_iterate
        end
    elseif Algorithm == :Newton
        if ((CompletedIters + SimpleStartIndex) % 2 == 1) & (CompletedIters > 1)
            # If we are in 3rd, 6th, 9th, 12th iterate from when we started Newton Acceleration then we want to do a Newton Iterate,
            # First we extract the guess that started this run of 3 iterates (x), the Function applied to it (fx) and the function applied to that (ffx)
            xk1  = Inputs[ :,(CompletedIters-1)]
            fxk1 = Outputs[:,(CompletedIters-1)]
            gxk1 = fxk1 .- xk1
            xk   = Inputs[ :,CompletedIters]
            fxk  = Outputs[:,CompletedIters]
            gxk  = fxk .- xk
            #ffx = Outputs[,CompletedIters ]
            # Now using the appropriate formula to make a new guess. Note that if a vector is input here it is used elementwise.
            derivative = (gxk.-gxk1)./(xk .- xk1)
            proposed_input   = xk .- (gxk./derivative)
        else
            # We just do a simple iterate.
            proposed_input = simple_iterate
        end
    elseif (Algorithm == :MPE) | (Algorithm == :RRE)
        SimpleIteratesMatrix = put_together_without_jumps(Inputs, Outputs)
        if (size(SimpleIteratesMatrix)[2] % ExtrapolationPeriod == 0)
            proposed_input = PolynomialExtrapolation(SimpleIteratesMatrix,Algorithm)
        else
            # We just do a simple iterate.
            proposed_input = simple_iterate
        end
    elseif (Algorithm == :VEA) | (Algorithm == :SEA)
        SimpleIteratesMatrix = put_together_without_jumps(Inputs, Outputs)
        if (size(SimpleIteratesMatrix)[2] % ExtrapolationPeriod) == 0
            proposed_input = EpsilonExtrapolation(SimpleIteratesMatrix, Algorithm)
        else
            # We just do a simple iterate.
            proposed_input = simple_iterate
        end
    else
        error("The algorithm you tried to input is not valid. Choose from :Simple, :Anderson, :Aitken, :Newton, :MPE, :RRE, :VEA or :SEA. Note capitalisation must match.")
    end
    # Now the replacement strategies - handle complex numbers properly
    if eltype(proposed_input) <: Complex
        dodgy_entries = (isnan.(real.(proposed_input)) .| isnan.(imag.(proposed_input))) .|
                       ismissing.(proposed_input) .|
                       (isinf.(real.(proposed_input)) .| isinf.(imag.(proposed_input)))
    else
        dodgy_entries = isnan.(proposed_input) .| ismissing.(proposed_input) .| isinf.(proposed_input)
    end

    if sum(dodgy_entries) != 0
        if ReplaceInvalids == :ReplaceElements
            proposed_input[dodgy_entries] = simple_iterate[dodgy_entries]
        elseif ReplaceInvalids == :ReplaceVector
            proposed_input = simple_iterate
        end
    end
    return @. (Dampening .* proposed_input) .+ ((1-Dampening) .* (Dampening_With_Input ? Inputs[:,end] : simple_iterate))
end
