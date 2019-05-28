"""
    execute_function_safely(Func::Function, x::Array{Float64,1})
This function creates a function that executes the function for which a fixed point is sought. It is a helper function that is not exported.
### Takes
 * Func - The function input to fixed_point
 * x - The point at which to evaluate the function.
### Returns
A FunctionEvaluationResult containing the following fields:
* Input_ - The input
* Output_ - The output of the FunctionExecutor. May be missing if function could not complete without error.
* Error_ - A enum of type FP_FunctionEvaluationError representing what error occured.
### Examples
Func(x) = sqrt(x)
execute_function_safely(Func, [-1.0,0.0,1.0])
execute_function_safely(Func,[Missing(),0.0,1.0])
execute_function_safely(Func,[7.0,0.0,1.0])
execute_function_safely(Func,[NaN,0.0,1.0])
execute_function_safely(Func,[Inf,0.0,1.0])
execute_function_safely(Func,-1.0)
execute_function_safely(Func,Missing())
execute_function_safely(Func,1.0)
execute_function_safely(Func,NaN)
execute_function_safely(Func,Inf)
"""
function execute_function_safely(Func::Function, x::Array{T,1}; type_check::Bool = false) where T<:Real
    # Check input
    if sum(ismissing.(x)) > 0
        return FunctionEvaluationResult(x, missing, InputMissingsDetected)
    elseif sum(isnan.(x)) > 0
        return FunctionEvaluationResult(x, missing, InputNAsDetected)
    elseif sum(isinf.(x)) > 0
        return FunctionEvaluationResult(x, missing, InputInfsDetected)
    end
    # Run function.
    lenx = length(x)
    tf_result = Array{T,1}(undef,lenx)
    try
        tf_result =  Func(x)
    catch
        return FunctionEvaluationResult(x, missing, ErrorExecutingFunction)
    end
    # Check Output and return.
    if sum(ismissing.(tf_result)) > 0
        return FunctionEvaluationResult(x, tf_result, OutputMissingsDetected)
    elseif sum(isnan.(tf_result)) > 0
        return FunctionEvaluationResult(x, tf_result, OutputNAsDetected)
    elseif sum(isinf.(tf_result)) > 0
        return FunctionEvaluationResult(x, tf_result, OutputInfsDetected)
    elseif (length(tf_result) != length(x))
        return FunctionEvaluationResult(x, tf_result, LengthOfOutputNotSameAsInput)
    elseif type_check & (!isa(tf_result[1], T))
        return FunctionEvaluationResult(x, tf_result, FunctionIsNotTypeStable)
    else
        return FunctionEvaluationResult(x, tf_result, NoError)
    end
end

"""
    fixed_point(func::Function, previous_FixedPointResults::FixedPointResults;
                    Algorithm::FixedPointAccelerationAlgorithm = Anderson,  ConvergenceMetric::Function  = supnorm(input::Array, output::Array) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::R = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::InvalidReplacement = NoAction, ConditionNumberThreshold::Real = 1e3) where R<:Real
    fixed_point(func::Function, Inputs::Array{T, 1};
                    Algorithm::FixedPointAccelerationAlgorithm = Anderson,  ConvergenceMetric::Function  = supnorm(input::Array, output::Array) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::InvalidReplacement = NoAction, ConditionNumberThreshold::Real = 1e3) where T<:Real
    fixed_point(func::Function, Inputs::Real;
                    Algorithm::FixedPointAccelerationAlgorithm = Anderson,  ConvergenceMetric::Function  = supnorm(input::Array, output::Array) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::InvalidReplacement = NoAction, ConditionNumberThreshold::Real = 1e3)
    fixed_point(func::Function, Inputs::Array{T, 2}; Outputs::Array{T,2} = Array{T,2}(undef,size(Inputs)[1],0),
                    Algorithm::FixedPointAccelerationAlgorithm = Anderson,  ConvergenceMetric::Function  = supnorm(input::Array, output::Array) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::InvalidReplacement = NoAction, ConditionNumberThreshold::Real = 1e3) where T<:Real where R<:Real

A function for finding the fixed point of another function
### Takes
 *  func - This is the function for which a fixed point is sought. This function must take and return a vector of the same size dimension.
 *  Inputs - This can be either a 1D-vector of values that is an initial guess for a fixed point or it can be an N x A matrix of previous inputs for which corresponding outputs are available. In this case N is the dimensionality of the fixed point vector you are seeking (Hence each column is a matrix that is input to f) and A is the number of previous Inputs/Outputs that are being provided to the fixed point. Where a matrix is input, a corresponding outputs must be provided or the last column of the outputs matrix is taken as a startpoint guess and the rest of the inputs and output matrices are discarded.
 *  Outputs - This is a matrix of the Function values for each column of the input. It must be provided so that column k of the outputs matrix is equal to Function(Column k of inputs matrix).
 *  Algorithm - This is the fixed point Algorithm to be used. It can be Anderson, Simple, Aitken, Newton, MPE, RRE, VEA or SEA. See vignette and references to see explanations of these Algorithms.
 *  ConvergenceMetric This is a function that takes in a vector of inputs and a table of outputs and returns a scalar. This scalar should be low when convergence is close to being achieved. By default this is the maximum residual by absolute value (the sup norm in the space of residuals).
 *   ConvergenceMetricThreshold This is the threshold for terminating the algorithm. The algorithm will terminate when the scalar that ConvergenceMetric returns is less than ConvergenceMetricThreshold. This can be set to a negative number in which case the algorithm will run until MaxIter is hit or an error occurs (Note that an error is likely in trying to use any Algorithm other than "Simple" when a fixed point is already found).
 *  MaxIter - This is the maximum number of iterates that will be undertaken.
 *  MaxM - This is the maximum number of saved iterates that are used in the Anderson algorithm. It has no effect if another Algorithm is chosen. Note that the number of previous iterates that will actually be used is the minimum of MaxIter, the dimensionality of the f's vector and the number of inputs that have been tried to  previously (the width of the Outputs matrix at each given stage of the algorithm). If PrintReports = TRUE, the number of previous iterates actually used is reported as the algorithm is running.
 *  ExtrapolationPeriod - This is the number of simple iterates to perform before extrapolating. This is used for the MPE, RRE, VEA and SEA Algorithms and has no effect if another Algorithm is chosen Where an epsilon algorithm is used this should be a multiple of 2, ie (4,6,8,etc).
 *  Dampening - This is the dampening parameter. By default it is 1 which means no dampening takes place. It can also be less than 1 (indicating dampening) or more than 1 (indicating extrapolation).
 *  PrintReports - This is a boolean describing whether to print ongoing ConvergenceMetric values for each iterate.
 *  ReportingSigFig - This is the number of significant figures that will be used in printing the convergence values to the console (only if PrintReports is TRUE).
 *  ReplaceInvalids -Sometimes an acceleration algorithm proposed a vector with an invalid coordinate (NaN, Inf or missing). This parameter can be set to ReplaceInvalids (to replace invalid coordinates by the simple iterate values), ReplaceVector (to replace entire vector with a simple iterate) or NoAction (where an imminent error will occur).
 *  ConditionNumberThreshold - This is a threshold for what condition number is acceptable for solving the least squares problem for the Anderson Algorithm. If the condition number is larger than this threshold then fewer previous iterates will be used in solving the problem. This has no effect unless the "Anderson" Algorithm is used.
### Returns
 * A list containing the fixed_point, the Inputs and corresponding Outputs, and convergence values (which are computed under the "ConvergenceMetric").
   The list will also include a "Finish" statement describing why it has finished. This is often going to be due to either MaxIter or ConvergenceMetricThreshold being
   reached. It may also terminate due to an error in generating a new input guess or using the function with that guess. If this occurs the function will terminate early
   and the "Finish" statement will describe the issue. In this event there will also be additional objects returned in the list "NewInputVector" and possibly
   "NewOutputVector" that are useful in debugging the issue.
### Examples
 #' # For the simplest possible example we can seek the fixed point of the cos function with a scalar.
 #' Inputs = 0.3
 #' Func(x) = cos(x)
 #' A = fixed_point(Func, Inputs; Algorithm = Aitken, Dampening = 0.5)
 #' B = fixed_point(Func, Inputs; Algorithm = Anderson, Dampening = 1.0)
 #'
 #' # For this next one the ConvergenceMetricThreshold is negative so the algorithm
 #' # will keep running until MaxIter is met.
 #' C = fixed_point(Func, Inputs; Algorithm = Simple, MaxIter = 4, ConvergenceMetricThreshold = -1.0)
 #' # But we can continue solving for this fixed point but now switching to the Newton Algorithm.
 #' D = fixed_point(Func, C[:Inputs], C[:Outputs]; Algorithm = Newton)
 #'
 #' # We can also find a 4 dimensional fixed point vector of this function.
 #' Inputs = [0.3, 98, 0, pi]
 #' E = fixed_point(Func, Inputs; Algorithm = Anderson)
 #' F = fixed_point(Func, Inputs; Algorithm = Anderson, MaxM = 4, ReportingSigFig = 13)
"""
function fixed_point(func::Function, previous_FixedPointResults::FixedPointResults;
                    Algorithm::FixedPointAccelerationAlgorithm = Anderson,  ConvergenceMetric::Function  = supnorm(input, output) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::R = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::InvalidReplacement = NoAction, ConditionNumberThreshold::Real = 1e3) where R<:Real
    Inputs = previous_FixedPointResults.Inputs_
    Outputs = previous_FixedPointResults.Outputs_
    return fixed_point(func, Inputs; Outputs = Outputs, Algorithm = Algorithm, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold,
                       MaxIter = MaxIter, MaxM = MaxM, ExtrapolationPeriod = ExtrapolationPeriod, Dampening = Dampening, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig,
                       ReplaceInvalids = ReplaceInvalids, ConditionNumberThreshold = ConditionNumberThreshold)
end
function fixed_point(func::Function, Inputs::Array{T, 1};
                    Algorithm::FixedPointAccelerationAlgorithm = Anderson,  ConvergenceMetric::Function  = supnorm(input, output) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::InvalidReplacement = NoAction, ConditionNumberThreshold::Real = 1e3) where T<:Real
    Inputs2 = Array{Float64, 2}(undef,size(Inputs)[1],1)
    Inputs2[:,1] = Inputs
    return fixed_point(func, Inputs2; Algorithm = Algorithm, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold,
                       MaxIter = MaxIter, MaxM = MaxM, ExtrapolationPeriod = ExtrapolationPeriod, Dampening = Dampening, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig,
                       ReplaceInvalids = ReplaceInvalids, ConditionNumberThreshold = ConditionNumberThreshold)
end
function fixed_point(func::Function, Inputs::Real;
                    Algorithm::FixedPointAccelerationAlgorithm = Anderson,  ConvergenceMetric::Function  = supnorm(input, output) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::InvalidReplacement = NoAction, ConditionNumberThreshold::Real = 1e3)
    Inputs2 = Array{typeof(Inputs), 2}(undef,1,1)
    Inputs2[1,1] = Inputs
    return fixed_point(func, Inputs2; Algorithm = Algorithm, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold,
                       MaxIter = MaxIter, MaxM = MaxM, ExtrapolationPeriod = ExtrapolationPeriod, Dampening = Dampening, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig,
                       ReplaceInvalids = ReplaceInvalids, ConditionNumberThreshold = ConditionNumberThreshold)
end
function fixed_point(func::Function, Inputs::Array{T, 2}; Outputs::Array{T,2} = Array{T,2}(undef,size(Inputs)[1],0),
                    Algorithm::FixedPointAccelerationAlgorithm = Anderson,  ConvergenceMetric::Function  = supnorm(input, output) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::InvalidReplacement = NoAction, ConditionNumberThreshold::Real = 1e3) where T<:Real where R<:Real
    # This code first tests if the input point is a fixed point. Then if it is not a while loop runs to try to find a fixed point.
    if (ConditionNumberThreshold < 1) error("ConditionNumberThreshold must be at least 1.")  end
    SimpleStartIndex = Integer(size(Outputs)[2])
    if isempty(Outputs)
        if size(Inputs)[2] > 1
            @warn("If you do not give outputs to the function then you can only give one vector of inputs (in a 2d array) to the fixed_pointFunction. So for a function that takes an N dimensional array you should input a Array{Float64}(N,1) array.  As you have input an array of size Array{Float64}(N,k) with k > 1 we have discarded everything but the last column to turn it into a Array{Float64}(N,1) array.\n")
            Inputs = Inputs[:,size(Inputs)[2]]
        end
    else
        if size(Inputs) != size(Outputs)
            @warn("If you input a matrix of outputs as well as a matrix of inputs then inputs and outputs must be the same shape. As they differ in this case the last column of the inputs matrix has been
                  taken as the starting point and everything else discarded.")
            Inputs  = Inputs[:,size(Inputs)[2]]
            Outputs = Array{T,2}(undef,size(Inputs)[1],0)
            SimpleStartIndex = Integer(1)
        end
    end
    LengthOfArray = size(Inputs)[1]
    output_type = T
    # Do an initial run if no runs have been done:
    if isempty(Outputs)
        ExecutedFunction = execute_function_safely(func, Inputs[:,1])
        if ExecutedFunction.Error_ != NoError
            return FixedPointResults(Inputs, Outputs, InvalidInputOrOutputOfIteration; FailedEvaluation_ = ExecutedFunction)
        end
        output_type = promote_type(typeof(Inputs[1]), typeof(ExecutedFunction.Output_[1]))
        converted_outputs = convert.(Ref(output_type), ExecutedFunction.Output_)
        Outputs = hcat(Outputs, converted_outputs)
        Inputs = convert.(Ref(output_type), Inputs)
    else
        # This ensures that MaxIter refers to max iter excluding any previous passed in results
        MaxIter = MaxIter + size(Outputs)[2]
        # This is to take into account previously passed in simple iterations (without jumps).
        SimpleStartIndex = SimpleStartIndex - (size(put_together_without_jumps(Inputs, Outputs))[2])
    end
    # First running through the last column of Inputs to test if we already have a fixed point.
    iter = Integer(size(Outputs)[2])
    ConvergenceVector = Array{output_type,1}(undef,iter)
    for i in 1:iter
        ConvergenceVector[i] = ConvergenceMetric(Inputs[:,i], Outputs[:,i])
    end
    if ConvergenceVector[iter] < ConvergenceMetricThreshold
        if (PrintReports)
            println("The last column of Inputs matrix is already a fixed point under input convergence metric and convergence threshold")
        end
        return FixedPointResults(Inputs, Outputs, ReachedConvergenceThreshold ; ConvergenceVector_  = vec(ConvergenceVector))
    end
    # Printing a report for initial convergence
    Convergence = ConvergenceVector[iter]
    if (PrintReports)
        println("                                          Algorithm: ", lpad(Algorithm, 8)   , ". Iteration: ", lpad(iter, 5),". Convergence: ", lpad(round(Convergence, digits=ReportingSigFig),ReportingSigFig+3))
    end
    iter = iter + 1
    while (Convergence > ConvergenceMetricThreshold) & (iter <= MaxIter)
        # Generating new input and output.
        NewInputFunctionReturn = fixed_point_new_input(Inputs, Outputs, Algorithm; MaxM = MaxM,
                                     SimpleStartIndex = SimpleStartIndex,ExtrapolationPeriod = ExtrapolationPeriod,
                                     Dampening = Dampening, ConditionNumberThreshold = ConditionNumberThreshold,
                                     PrintReports = PrintReports, ReplaceInvalids = ReplaceInvalids)
        if PrintReports & (Algorithm != Anderson)
            print(lpad("",42))
        end
        ExecutedFunction = execute_function_safely(func, NewInputFunctionReturn; type_check = true)
        if ExecutedFunction.Error_ != NoError
            return FixedPointResults(Inputs, Outputs, InvalidInputOrOutputOfIteration; ConvergenceVector_  = vec(ConvergenceVector), FailedEvaluation_ = ExecutedFunction)
        end
        Inputs  = hcat(Inputs, ExecutedFunction.Input_)
        Outputs = hcat(Outputs, convert(Array{output_type,1}, ExecutedFunction.Output_))
        # Checking and recording convergence
        Convergence = ConvergenceMetric(ExecutedFunction.Input_, ExecutedFunction.Output_)
        ConvergenceVector =  vcat(ConvergenceVector, Convergence)
        # Output of report and going to next iteration.
        if (PrintReports) println("Algorithm: ", lpad(Algorithm,8)   , ". Iteration: ", lpad(iter,5), ". Convergence: ", lpad(round(Convergence, digits=ReportingSigFig),ReportingSigFig+3)) end
        iter  = iter + 1
    end
    Finish = ReachedMaxIter
    if (Convergence < ConvergenceMetricThreshold) Finish = ReachedConvergenceThreshold end
    return FixedPointResults(Inputs, Outputs, Finish; ConvergenceVector_  = vec(ConvergenceVector))
end

"""
    fixed_point_new_input(Inputs::AbstractArray{T,2}, Outputs::AbstractArray{T,2}, Algorithm::FixedPointAccelerationAlgorithm = Anderson;
                               MaxM::Integer = 10, SimpleStartIndex::Integer = 1, ExtrapolationPeriod::Integer = 1, Dampening::S = AbstractFloat(1),
                               ConditionNumberThreshold::R = AbstractFloat(1000), PrintReports::Bool = false, ReplaceInvalids::InvalidReplacement = NoAction) where R<:Real where S<:Real where T<:Real

This function takes the previous inputs and outputs from the fixed_point function and determines what vector to try next in seeking a fixed point.
### Takes
 *  Inputs - This is an N x A matrix of previous inputs for which corresponding outputs are available. In this case N is the dimensionality of the fixed point vector that is being sought (Hence each column is a matrix that is input to the "Function") and A is the number of previous Inputs/Outputs that are being provided to the fixed point.
* Outputs - This is a matrix of Function values for the each column of the "Inputs" matrix.
* Algorithm - This is the fixed point Algorithm to be used. It can be "Anderson", "Simple", "Aitken", "Newton", "MPE", "RRE", "VEA", "SEA".
* MaxM - This is the number of saved iterates that are used in the Anderson algorithm. This has no role if another Algorithm is used.
* SimpleStartIndex - This is the index for what column in the input/output matrices did the algorithm start doing simple iterates without jumps. This is used for all Algorithms except the simple and Anderson Algorithms where it has no effect.
* ExtrapolationPeriod - This is the number of simple iterates to perform before extrapolating. This is used for the MPE, RRE, VEA and SEA Algorithms and has no effect if another Algorithm is chosen.
* Dampening - This is the dampening parameter. By default it is 1 which means no dampening takes place. It can also be less than 1 (indicating dampening) or more than 1 (indicating extrapolation).
* ConditionNumberThreshold - This is what threshold should be chosen to drop previous iterates if the matrix is ill conditioned. Only used in Anderson acceleration.
* PrintReports - This is a boolean describing whether to print ongoing ConvergenceMetric values for each iterate.
### Returns
 * A nicely formatted string version of the input number for printing to the console.
### Examples
FPFunction = function(x){c(0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2])}
A = fixed_point( Function = FPFunction, Inputs = [0.3,900], MaxIter = 6, Algorithm = Simple)
NewGuessAnderson = fixed_point_new_input(A[:Inputs], A[:Outputs], Algorithm = Anderson)
NewGuessVEA = fixed_point_new_input(A[:Inputs], A[:Outputs], Algorithm = VEA)
NewGuessMPE = fixed_point_new_input(A[:Inputs], A[:Outputs], Algorithm = MPE)
NewGuessAitken = fixed_point_new_input(A[:Inputs], A[:Outputs], Algorithm = Aitken)
"""
function fixed_point_new_input(Inputs::AbstractArray{T,2}, Outputs::AbstractArray{T,2}, Algorithm::FixedPointAccelerationAlgorithm = Anderson;
                               MaxM::Integer = 10, SimpleStartIndex::Integer = 1, ExtrapolationPeriod::Integer = 1, Dampening::S = AbstractFloat(1),
                               ConditionNumberThreshold::R = AbstractFloat(1000), PrintReports::Bool = false, ReplaceInvalids::InvalidReplacement = NoAction) where R<:Real where S<:Real where T<:Real
    CompletedIters = size(Outputs)[2]
    simple_iterate = Outputs[:,CompletedIters]
    proposed_input = repeat([NaN], size(simple_iterate)[1])
    if Algorithm == Simple
         proposed_input = simple_iterate
    elseif Algorithm == Anderson
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
            while sum(isnan.(Coeffs)) > 0
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
                Fit = fit(LinearModel,  hcat(DeltaResids), LastResid)
                Coeffs = Fit.pp.beta0
                if sum(isnan.(Coeffs)) > 0
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
                if (PrintReports) print("Condition number is ", lpad(ConditionNumber, 5),". Used:",  lpad(M+1, 3)," lags. ") end
                proposed_input = LastOutput .- (Dampening .* vec(DeltaOutputs * Coeffs))
            end
        end
    elseif Algorithm == Aitken
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
    elseif Algorithm == Newton
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
    elseif (Algorithm == MPE) | (Algorithm == RRE)
        SimpleIteratesMatrix = put_together_without_jumps(Inputs, Outputs)
        if (size(SimpleIteratesMatrix)[2] % ExtrapolationPeriod == 0)
            proposed_input = PolynomialExtrapolation(SimpleIteratesMatrix,Algorithm)
        else
            # We just do a simple iterate.
            proposed_input = simple_iterate
        end
    elseif (Algorithm == VEA) | (Algorithm == SEA)
        SimpleIteratesMatrix = put_together_without_jumps(Inputs, Outputs)
        if (size(SimpleIteratesMatrix)[2] % ExtrapolationPeriod) == 0
            proposed_input = EpsilonExtrapolation(SimpleIteratesMatrix, Algorithm)
        else
            # We just do a simple iterate.
            proposed_input = simple_iterate
        end
    end
    # Now the replacement strategies
    dodgy_entries = isnan.(proposed_input) .| ismissing.(proposed_input) .| isinf.(proposed_input)
    if sum(dodgy_entries) != 0
        if ReplaceInvalids == ReplaceElements
            proposed_input[dodgy_entries] = simple_iterate[dodgy_entries]
        elseif ReplaceInvalids == ReplaceVector
            proposed_input = simple_iterate
        end
    end
    return @. (Dampening .* proposed_input) .+ ((1-Dampening) .* simple_iterate)
end

"""
This function performs Minimal Polynomial extrapolation (MPE) or Reduced Rank Extrapolation (RRE) given a matrix of previous iterates of the function.
### Takes
 * Iterates - A matrix of inputs and outputs excluding jumps. Can be pieced together from Inputs and Outputs matrices of the fixed_point function using the put_together_without_jumps function
 * Algorithm The Algorithm for polynomial extrapolation. Should be either "MPE" for minimal polynomial extrapolation or "RRE" for reduced rank extrapolation.
### Returns
 * A vector containing the extrapolated vector.
"""
function PolynomialExtrapolation(Iterates::AbstractArray{R,2}, Algorithm::FixedPointAccelerationAlgorithm) where R<:Real
    if (Algorithm == MPE)
        TotalColumnsOfIterates = size(Iterates)[2]
        OldDifferences         = Iterates[:,2:(TotalColumnsOfIterates-1)] .- Iterates[:,1:(TotalColumnsOfIterates-2)]
        LastDifference         = Iterates[:,TotalColumnsOfIterates]       .- Iterates[:,(TotalColumnsOfIterates-1)]
        InverseOldDifferences  = pinv(OldDifferences)
        cVector                = -InverseOldDifferences * LastDifference
        cVector                = vcat(cVector,1)
        sumVec                 = sum(cVector)
        return (Iterates[:,2:TotalColumnsOfIterates] * cVector) ./ sumVec
    elseif (Algorithm == RRE)
        TotalColumnsOfIterates   = size(Iterates)[2]
        FirstColumn              = Iterates[:,1]
        Differences              = Iterates[:,2:(TotalColumnsOfIterates)]      - Iterates[:,1:(TotalColumnsOfIterates-1)]
        SecondDifferences        = Differences[:,2:(TotalColumnsOfIterates-1)] - Differences[:,1:(TotalColumnsOfIterates-2)]
        FirstDifference          = Differences[:,1]
        Differences              = Differences[:,1:(TotalColumnsOfIterates-2)]
        InverseSecondDifferences = pinv(SecondDifferences)
        return FirstColumn - ((Differences * InverseSecondDifferences) * FirstDifference)
    else
        error("Invalid Algorithm input. PolynomialExtrapolation function can only take Algorithm as MPE or RRE.")
    end
end

"""
This is a helper function for EpsilonExtrapolation
### Takes
 * Iterates - A matrix representing different iterates with one iterate per column. Can be pieced together from Inputs and Outputs matrices of the fixed_point function using the put_together_without_jumps function
 * Algorithm - Algorithm for epsilon extrapolation. Should be either "VEA" for the vector extrapolation algorithm or "SEA" for the scalar epsilon algorithm.
### Returns
 *  A vector with the extrapolated vector.
"""
function EpsilonExtrapolation(Iterates::AbstractArray{R,2}, Algorithm::FixedPointAccelerationAlgorithm) where R<:Real
    # The function cannot do anything to a one column input so will return input unchanged.
    if (size(Iterates)[2] == 1) return Iterates end
    if (size(Iterates)[2] % 2 == 0) Iterates = Iterates[:,2:size(Iterates)[2]] end
    if (!(Algorithm in [VEA, SEA])) error("Invalid Algorithm input. EpsilonExtrapolation function can only take Algorithm as VEA or SEA") end
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
    if (any(isnan.(Mat)) | any(ismissing.(Mat)))
        Mat = EpsilonExtrapolation(Iterates[:,3:size(Iterates)[2]],Algorithm)
    end
    return Mat[:,1]
end

"""
This is a helper function for EpsilonExtrapolation
### Takes
 * DifferenceMatrix - The matrix of the differences in elements to be inverted.
 * Algorithm - SEA or VEA.
### Returns
 * A vector of the result of inverting each (column) vector in a mmatrix.
"""
function EpsilonExtrapolationVectorOfInverses(DifferenceMatrix::AbstractArray{T,2}, Algorithm::FixedPointAccelerationAlgorithm) where T<:Real
    if (size(DifferenceMatrix)[1] == 1) | (Algorithm == SEA)
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

"""
This function takes the previous inputs and outputs and assembles a matrix with both excluding jumps.
### Takes
 * Inputs - This is an N x A matrix of previous inputs for which corresponding outputs are available. In this case N is the dimensionality of the fixed point vector that is being sought (and each column is a matrix that is input to the "Function") and A is the number of previous Inputs/Outputs that are being provided to the fixed point.
 * Outputs This is a matrix of "Function" values for each column of the "Inputs" matrix.
 * AgreementThreshold A parameter for determining when a column in Inputs and a column in Outputs match. They are deemed to match if the sum of the absolute values of the difference in the columns is less than AgreementThreshold.
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
