"""
    fixed_point(func::Function, previous_FixedPointResults::FixedPointResults;
                    Algorithm::Symbol = :Anderson,  ConvergenceMetric::Function  = supnorm(input::Array, output::Array) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::R = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::Symbol = :NoAction, ConditionNumberThreshold::Real = 1e3, quiet_errors::Bool = false) where R<:Real
    fixed_point(func::Function, Inputs::Array{T, 1};
                    Algorithm::Symbol = :Anderson,  ConvergenceMetric::Function  = supnorm(input::Array, output::Array) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::Symbol = :NoAction, ConditionNumberThreshold::Real = 1e3, quiet_errors::Bool = false) where T<:Real
    fixed_point(func::Function, Inputs::Real;
                    Algorithm::Symbol = :Anderson,  ConvergenceMetric::Function  = supnorm(input::Array, output::Array) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::Symbol = :NoAction, ConditionNumberThreshold::Real = 1e3, quiet_errors::Bool = false)
    fixed_point(func::Function, Inputs::Array{T, 2}; Outputs::Array{T,2} = Array{T,2}(undef,size(Inputs)[1],0),
                    Algorithm::Symbol = :Anderson,  ConvergenceMetric::Function  = supnorm(input::Array, output::Array) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Real = AbstractFloat(1.0),
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::Symbol = :NoAction, ConditionNumberThreshold::Real = 1e3, quiet_errors::Bool = false, other_outputs::Union{Missing,NamedTuple} = missing) where T<:Real where R<:Real

A function for finding the fixed point of another function
### Inputs
 *  `func` - This is the function for which a fixed point is sought. This function must take and return a vector of the same size dimension.
 *  `Inputs` - This can be either a 1D-vector of values that is an initial guess for a fixed point or it can be an N x A matrix of previous inputs for which corresponding outputs are available. In this case N is the dimensionality of the fixed point vector you are seeking (Hence each column is a matrix that is input to f) and A is the number of previous Inputs/Outputs that are being provided to the fixed point. Where a matrix is input, a corresponding outputs must be provided or the last column of the outputs matrix is taken as a startpoint guess and the rest of the inputs and output matrices are discarded.
 *  `Outputs` - This is a matrix of the Function values for each column of the input. It must be provided so that column k of the outputs matrix is equal to Function(Column k of inputs matrix).
 *  `Algorithm` - This is the fixed point Algorithm to be used. It can be :Anderson, :Simple, :Aitken, :Newton, :MPE, :RRE, :VEA or :SEA. See documentation and references to see explanations of these Algorithms.
 *  `ConvergenceMetric` - This is a function that takes in a vector of inputs and a table of outputs and returns a scalar. This scalar should be low when convergence is close to being achieved. By default this is the maximum residual by absolute value (the sup norm in the space of residuals).
 *  `ConvergenceMetricThreshold` - This is the threshold for terminating the algorithm. The algorithm will terminate when the scalar that ConvergenceMetric returns is less than ConvergenceMetricThreshold. This can be set to a negative number in which case the algorithm will run until MaxIter is hit or an error occurs (Note that an error is likely in trying to use any Algorithm other than "Simple" when a fixed point is already found).
 *  `MaxIter` - This is the maximum number of iterates that will be undertaken.
 *  `MaxM` - This is the maximum number of saved iterates that are used in the Anderson algorithm. It has no effect if another Algorithm is chosen. Note that the number of previous iterates that will actually be used is the minimum of MaxIter, the dimensionality of the f's vector and the number of inputs that have been tried to  previously (the width of the Outputs matrix at each given stage of the algorithm). If PrintReports = TRUE, the number of previous iterates actually used is reported as the algorithm is running.
 *  `ExtrapolationPeriod` - This is the number of simple iterates to perform before extrapolating. This is used for the MPE, RRE, VEA and SEA Algorithms and has no effect if another Algorithm is chosen Where an epsilon algorithm is used this should be a multiple of 2, ie (4,6,8,etc).
 *  `Dampening` - This is the dampening parameter. By default it is 1 which means no dampening takes place. It can also be less than 1 (indicating dampening) or more than 1 (indicating extrapolation).
 * `Dampening_With_Input` - This is a boolean that indicates whether the dampening parameter should be multiplied by the input (if true) or the output of the most recent iterate.
 *  `PrintReports` - This is a boolean describing whether to print ongoing ConvergenceMetric values for each iterate.
 *  `ReportingSigFig` - This is the number of significant figures that will be used in printing the convergence values to the console (only if PrintReports is TRUE).
 *  `ReplaceInvalids` - Sometimes an acceleration algorithm proposed a vector with an invalid coordinate (NaN, Inf or missing). This parameter can be set to :ReplaceInvalids (to replace invalid coordinates by the simple iterate values), :ReplaceVector (to replace entire vector with a simple iterate) or :NoAction (where an imminent error will occur).
 *  `ConditionNumberThreshold` - This is a threshold for what condition number is acceptable for solving the least squares problem for the Anderson Algorithm. If the condition number is larger than this threshold then fewer previous iterates will be used in solving the problem. This has no effect unless the :Anderson Algorithm is used.
 *  `quiet_errors` - If true the function will return everything already calculated as soon as an error occurs. The callstack that lead to the error is not returned however. If false an error will be thrown with a callstack.
 *  `other_outputs` - This allows you to pass in side products (as in Other_Output_ in a FixedPointResults struct). It is only used if the FixedPointResults that is input has already found a fixedpoint.
### Returns
 * A `FixedPointResults` struct containing the fixed_point, the Inputs and corresponding Outputs, and convergence values (which are computed under the "ConvergenceMetric").
   The list will also include a "Finish" statement describing why it has finished. This is often going to be due to either MaxIter or ConvergenceMetricThreshold being
   reached. It may also terminate due to an error in generating a new input guess or using the function with that guess. If this occurs the function will terminate early
   and the "Finish" statement will describe the issue. In this event there will also be additional objects returned in the list "NewInputVector" and possibly
   "NewOutputVector" that are useful in debugging the issue.
### Examples
    # For the simplest possible example we can seek the fixed point of the cos function with a scalar.
    Inputs = 0.3
    Func(x) = cos(x)
    A = fixed_point(Func, Inputs; Algorithm = :Aitken, Dampening = 0.5)
    B = fixed_point(Func, Inputs; Algorithm = :Anderson, Dampening = 1.0)

    # For this next one the ConvergenceMetricThreshold is negative so the algorithm
    # will keep running until MaxIter is met.
    C = fixed_point(Func, Inputs; Algorithm = :Simple, MaxIter = 4, ConvergenceMetricThreshold = -1.0)
    # But we can continue solving for this fixed point but now switching to the Newton Algorithm.
    D = fixed_point(Func, C[:Inputs], C[:Outputs]; Algorithm = :Newton)

    # We can also find a 4 dimensional fixed point vector of this function.
    Inputs = [0.3, 98, 0, pi]
    E = fixed_point(Func, Inputs; Algorithm = :Anderson)
    F = fixed_point(Func, Inputs; Algorithm = :Anderson, MaxM = 4, ReportingSigFig = 13)
"""
function fixed_point(func::Function, previous_FixedPointResults::FixedPointResults;
                    Algorithm::Symbol = :Anderson,  ConvergenceMetric::Function  = supnorm(input, output) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Number = 1.0, Dampening_With_Input::Bool = false,
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::Symbol = :NoAction, ConditionNumberThreshold::Real = 1e3, quiet_errors::Bool = false)
    Inputs = previous_FixedPointResults.Inputs_
    Outputs = previous_FixedPointResults.Outputs_
    side_products = previous_FixedPointResults.Other_Output_
    return fixed_point(func, Inputs; Outputs = Outputs, Algorithm = Algorithm, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold,
                       MaxIter = MaxIter, MaxM = MaxM, ExtrapolationPeriod = ExtrapolationPeriod, Dampening = Dampening, Dampening_With_Input = Dampening_With_Input, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig,
                       ReplaceInvalids = ReplaceInvalids, ConditionNumberThreshold = ConditionNumberThreshold, quiet_errors = quiet_errors, other_outputs = side_products)
end

function fixed_point(func::Function, Inputs::Array{T, 1};
                    Algorithm::Symbol = :Anderson,  ConvergenceMetric::Function  = supnorm(input, output) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Number = 1.0, Dampening_With_Input::Bool = false,
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::Symbol = :NoAction, ConditionNumberThreshold::Real = 1e3, quiet_errors::Bool = false) where T<:Number
    Inputs2 = Array{T, 2}(undef,size(Inputs)[1],1)
    Inputs2[:,1] = Inputs
    return fixed_point(func, Inputs2; Algorithm = Algorithm, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold,
                       MaxIter = MaxIter, MaxM = MaxM, ExtrapolationPeriod = ExtrapolationPeriod, Dampening = Dampening, Dampening_With_Input = Dampening_With_Input, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig,
                       ReplaceInvalids = ReplaceInvalids, ConditionNumberThreshold = ConditionNumberThreshold, quiet_errors = quiet_errors)
end

function fixed_point(func::Function, Inputs::Number;
                    Algorithm::Symbol = :Anderson,  ConvergenceMetric::Function  = supnorm(input, output) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Number = 1.0, Dampening_With_Input::Bool = false,
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::Symbol = :NoAction, ConditionNumberThreshold::Real = 1e3, quiet_errors::Bool = false)
    Inputs2 = Array{typeof(Inputs), 2}(undef,1,1)
    Inputs2[1,1] = Inputs
    return fixed_point(func, Inputs2; Algorithm = Algorithm, ConvergenceMetric = ConvergenceMetric, ConvergenceMetricThreshold = ConvergenceMetricThreshold,
                       MaxIter = MaxIter, MaxM = MaxM, ExtrapolationPeriod = ExtrapolationPeriod, Dampening = Dampening, Dampening_With_Input = Dampening_With_Input, PrintReports = PrintReports, ReportingSigFig = ReportingSigFig,
                       ReplaceInvalids = ReplaceInvalids, ConditionNumberThreshold = ConditionNumberThreshold, quiet_errors = quiet_errors)
end

function fixed_point(func::Function, Inputs::Array{T, 2}; Outputs::Array{<:Number,2} = Array{T,2}(undef,size(Inputs)[1],0),
                    Algorithm::Symbol = :Anderson,  ConvergenceMetric::Function  = supnorm(input, output) = maximum(abs.(output .- input)),
                    ConvergenceMetricThreshold::Real = 1e-10, MaxIter::Integer = Integer(1000), MaxM::Integer = Integer(10), ExtrapolationPeriod::Integer = Integer(7), Dampening::Number = 1.0, Dampening_With_Input::Bool = false,
                    PrintReports::Bool = false, ReportingSigFig::Integer = Integer(10), ReplaceInvalids::Symbol = :NoAction, ConditionNumberThreshold::Real = 1e3, quiet_errors::Bool = false, other_outputs::Union{Missing,NamedTuple} = missing) where T<:Number
    # This code first tests if the input point is a fixed point. Then if it is not a while loop runs to try to find a fixed point.
    if (ConditionNumberThreshold < 1) error("ConditionNumberThreshold must be at least 1.")  end
    SimpleStartIndex = Integer(size(Outputs)[2])
    if isempty(Outputs)
        if size(Inputs)[2] > 1
            @warn("If you do not give outputs to the function then you can only give one vector of inputs (in a 2d array) to the fixed_pointFunction. So for a function that takes an N dimensional array you should input a Array{Float64}(N,1) array.  As you have input an array of size Array{Float64}(N,k) with k > 1 we have discarded everything but the last column to turn it into a Array{Float64}(N,1) array.\n")
            Inputs = Inputs[:,size(Inputs)[2]]
            Inputs = reshape(Inputs, length(Inputs), 1)
        end
    else
        if size(Inputs) != size(Outputs)
            @warn("If you input a matrix of outputs as well as a matrix of inputs then inputs and outputs must be the same shape. As they differ in this case the last column of the inputs matrix has been taken as the starting point and everything else discarded.")
            Inputs  = Inputs[:,size(Inputs)[2]]
            Inputs = reshape(Inputs, length(Inputs), 1)
            Outputs = Array{T,2}(undef,size(Inputs)[1],0)
            SimpleStartIndex = Integer(1)
        end
    end
    LengthOfArray = size(Inputs)[1]
    output_type = promote_type(T,eltype(Outputs))
    final_other_output = other_outputs
    # Do an initial run if no runs have been done:
    if isempty(Outputs)
        ExecutedFunction = execute_function_safely(func, Inputs[:,1]; quiet_errors = quiet_errors)
        final_other_output = ExecutedFunction.Other_Output_
        if ExecutedFunction.Error_ != :NoError
            return FixedPointResults(Inputs, Outputs, :InvalidInputOrOutputOfIteration; FailedEvaluation_ = ExecutedFunction, Other_Output = final_other_output)
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
    # Convergence metrics should always return real numbers, even for complex inputs
    convergence_type = output_type <: Complex ? real(output_type) : output_type
    ConvergenceVector = Array{convergence_type,1}(undef,iter)
    for i in 1:iter
        ConvergenceVector[i] = ConvergenceMetric(Inputs[:,i], Outputs[:,i])
    end
    if ConvergenceVector[iter] < ConvergenceMetricThreshold
        if (PrintReports)
            println("The last column of Inputs matrix is already a fixed point under input convergence metric and convergence threshold")
        end
        return FixedPointResults(Inputs, Outputs, :ReachedConvergenceThreshold; ConvergenceVector_ = vec(ConvergenceVector), Other_Output = final_other_output)
    end
    # Printing a report for initial convergence
    Convergence = ConvergenceVector[iter]
    if (PrintReports)
        println("                                          Algorithm: ", lpad(Algorithm, 8)   , ". Iteration: ", lpad(iter, 5),". Convergence: ", lpad(round(Convergence, sigdigits=ReportingSigFig),ReportingSigFig+4), ". Time: ", now())
    end
    iter = iter + 1
    while (Convergence > ConvergenceMetricThreshold) & (iter <= MaxIter)
        # Generating new input and output.
        NewInputFunctionReturn = fixed_point_new_input(Inputs, Outputs, Algorithm; MaxM = MaxM,
                                     SimpleStartIndex = SimpleStartIndex,ExtrapolationPeriod = ExtrapolationPeriod,
                                     Dampening = Dampening, Dampening_With_Input = Dampening_With_Input, ConditionNumberThreshold = ConditionNumberThreshold,
                                     PrintReports = PrintReports, ReplaceInvalids = ReplaceInvalids)
        if PrintReports & (Algorithm != :Anderson)
            print(lpad("",42))
        end
        ExecutedFunction = execute_function_safely(func, NewInputFunctionReturn; type_check = true, quiet_errors = quiet_errors)
        if ExecutedFunction.Error_ != :NoError
            return FixedPointResults(Inputs, Outputs, :InvalidInputOrOutputOfIteration; ConvergenceVector_  = vec(ConvergenceVector), FailedEvaluation_ = ExecutedFunction, Other_Output = ExecutedFunction.Other_Output_)
        end
        final_other_output = ExecutedFunction.Other_Output_
        Inputs  = hcat(Inputs, ExecutedFunction.Input_)
        Outputs = hcat(Outputs, convert(Array{output_type,1}, ExecutedFunction.Output_))
        # Checking and recording convergence
        Convergence = ConvergenceMetric(ExecutedFunction.Input_, ExecutedFunction.Output_)
        ConvergenceVector =  vcat(ConvergenceVector, Convergence)
        # Output of report and going to next iteration.
        if (PrintReports) println("Algorithm: ", lpad(Algorithm,8)   , ". Iteration: ", lpad(iter,5), ". Convergence: ", lpad(round(Convergence, sigdigits=ReportingSigFig),ReportingSigFig+4), ". Time: ", now()) end
        iter  = iter + 1
    end
    Finish = (Convergence < ConvergenceMetricThreshold) ? :ReachedConvergenceThreshold  : :ReachedMaxIter
    return FixedPointResults(Inputs, Outputs, Finish; ConvergenceVector_  = vec(ConvergenceVector),  Other_Output = final_other_output)
end
