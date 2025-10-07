"""
    fixed_point(func::Function, initial_guess, algorithm::FixedPointAlgorithm = Anderson();
                ConvergenceMetric::Function = (input, output) -> maximum(abs.(output .- input)),
                ConvergenceMetricThreshold::Real = 1e-10,
                MaxIter::Integer = 1000,
                Dampening::Real = 1.0,
                Dampening_With_Input::Bool = false,
                PrintReports::Bool = false,
                ReportingSigFig::Integer = 10,
                ReplaceInvalids::Symbol = :NoAction,
                quiet_errors::Bool = false)

    fixed_point(func::Function, previous_results::FixedPointResults, algorithm::FixedPointAlgorithm = Anderson();
                ConvergenceMetric::Function = (input, output) -> maximum(abs.(output .- input)),
                ConvergenceMetricThreshold::Real = 1e-10,
                MaxIter::Integer = 1000,
                Dampening::Real = 1.0,
                Dampening_With_Input::Bool = false,
                PrintReports::Bool = false,
                ReportingSigFig::Integer = 10,
                ReplaceInvalids::Symbol = :NoAction,
                quiet_errors::Bool = false)

    fixed_point(func::Function, inputs::Matrix, outputs::Matrix, algorithm::FixedPointAlgorithm = Anderson();
                ConvergenceMetric::Function = (input, output) -> maximum(abs.(output .- input)),
                ConvergenceMetricThreshold::Real = 1e-10,
                MaxIter::Integer = 1000,
                Dampening::Real = 1.0,
                Dampening_With_Input::Bool = false,
                PrintReports::Bool = false,
                ReportingSigFig::Integer = 10,
                ReplaceInvalids::Symbol = :NoAction,
                quiet_errors::Bool = false,
                other_outputs::Union{Missing,NamedTuple} = missing)

Find the fixed point of a function using various acceleration algorithms.

### Inputs
 *  `func` - The function for which a fixed point is sought. Must take and return a vector of the same size.
 *  `initial_guess` - Initial guess for the fixed point. Can be a scalar, vector, or matrix.
 *  `algorithm` - Algorithm instance (e.g., `Anderson()`, `Simple()`, `MPE(extrapolation_period=4)`).
 *  `previous_results` - Continue from previous `FixedPointResults`.
 *  `inputs` - N×A matrix of previous inputs (advanced usage).
 *  `outputs` - N×A matrix of corresponding outputs (advanced usage).
 *  `ConvergenceMetric` - Function that measures convergence (input, output) -> scalar.
 *  `ConvergenceMetricThreshold` - Threshold for convergence (default: 1e-10).
 *  `MaxIter` - Maximum number of iterations (default: 1000).
 *  `Dampening` - Dampening parameter (default: 1.0, no dampening).
 *  `Dampening_With_Input` - Apply dampening to input vs output (default: false).
 *  `PrintReports` - Print iteration progress (default: false).
 *  `ReportingSigFig` - Significant figures for progress reports (default: 10).
 *  `ReplaceInvalids` - How to handle NaN/Inf: `:NoAction`, `:ReplaceElements`, `:ReplaceVector`.
 *  `quiet_errors` - Return partial results on error instead of throwing (default: false).
 *  `other_outputs` - Additional outputs to pass through (advanced usage).
### Returns
 * A `FixedPointResults` struct containing the fixed_point, the Inputs and corresponding Outputs, and convergence values (which are computed under the "ConvergenceMetric").
   The list will also include a "Finish" statement describing why it has finished. This is often going to be due to either MaxIter or ConvergenceMetricThreshold being
   reached. It may also terminate due to an error in generating a new input guess or using the function with that guess. If this occurs the function will terminate early
   and the "Finish" statement will describe the issue. In this event there will also be additional objects returned in the list "NewInputVector" and possibly
   "NewOutputVector" that are useful in debugging the issue.
### Examples
    # Simple scalar function
    f(x) = cos(x)
    result = fixed_point(f, 0.3, Aitken(); Dampening = 0.5)

    # Vector function with Anderson acceleration
    g(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
    result = fixed_point(g, [0.3, 900.0], Anderson(maxM=5))

    # MPE with custom extrapolation period
    result = fixed_point(g, [0.3, 900.0], MPE(extrapolation_period=4))

    # Continue from previous results with different algorithm
    result2 = fixed_point(f, result, Simple())

    # Run for fixed iterations regardless of convergence
    result3 = fixed_point(g, [0.3, 900.0], Simple(); MaxIter=4, ConvergenceMetricThreshold=-1.0)
"""
# Main interface: continue from previous results
function fixed_point(
    func::Function,
    previous_results::FixedPointResults,
    algorithm::FixedPointAlgorithm = Anderson();
    ConvergenceMetric::Function = (input, output) -> maximum(abs.(output .- input)),
    ConvergenceMetricThreshold::Real = 1e-10,
    MaxIter::Integer = 1000,
    Dampening::Real = 1.0,
    Dampening_With_Input::Bool = false,
    PrintReports::Bool = false,
    ReportingSigFig::Integer = 10,
    ReplaceInvalids::Symbol = :NoAction,
    quiet_errors::Bool = false,
)
    return fixed_point(
        func,
        previous_results.Inputs_,
        previous_results.Outputs_,
        algorithm;
        ConvergenceMetric=ConvergenceMetric,
        ConvergenceMetricThreshold=ConvergenceMetricThreshold,
        MaxIter=MaxIter,
        Dampening=Dampening,
        Dampening_With_Input=Dampening_With_Input,
        PrintReports=PrintReports,
        ReportingSigFig=ReportingSigFig,
        ReplaceInvalids=ReplaceInvalids,
        quiet_errors=quiet_errors,
        other_outputs=previous_results.Other_Output_,
    )
end

# Main interface: vector input
function fixed_point(
    func::Function,
    initial_guess::AbstractVector{T},
    algorithm::FixedPointAlgorithm = Anderson();
    ConvergenceMetric::Function = (input, output) -> maximum(abs.(output .- input)),
    ConvergenceMetricThreshold::Real = 1e-10,
    MaxIter::Integer = 1000,
    Dampening::Real = 1.0,
    Dampening_With_Input::Bool = false,
    PrintReports::Bool = false,
    ReportingSigFig::Integer = 10,
    ReplaceInvalids::Symbol = :NoAction,
    quiet_errors::Bool = false,
) where {T<:Number}
    inputs_matrix = reshape(initial_guess, length(initial_guess), 1)
    outputs_matrix = Array{T,2}(undef, size(inputs_matrix)[1], 0)
    return fixed_point(
        func,
        inputs_matrix,
        outputs_matrix,
        algorithm;
        ConvergenceMetric=ConvergenceMetric,
        ConvergenceMetricThreshold=ConvergenceMetricThreshold,
        MaxIter=MaxIter,
        Dampening=Dampening,
        Dampening_With_Input=Dampening_With_Input,
        PrintReports=PrintReports,
        ReportingSigFig=ReportingSigFig,
        ReplaceInvalids=ReplaceInvalids,
        quiet_errors=quiet_errors,
    )
end

# Main interface: scalar input
function fixed_point(
    func::Function,
    initial_guess::Number,
    algorithm::FixedPointAlgorithm = Anderson();
    ConvergenceMetric::Function = (input, output) -> maximum(abs.(output .- input)),
    ConvergenceMetricThreshold::Real = 1e-10,
    MaxIter::Integer = 1000,
    Dampening::Real = 1.0,
    Dampening_With_Input::Bool = false,
    PrintReports::Bool = false,
    ReportingSigFig::Integer = 10,
    ReplaceInvalids::Symbol = :NoAction,
    quiet_errors::Bool = false,
)
    inputs_matrix = Array{typeof(initial_guess),2}(undef, 1, 1)
    inputs_matrix[1, 1] = initial_guess
    outputs_matrix = Array{typeof(initial_guess),2}(undef, 1, 0)
    return fixed_point(
        func,
        inputs_matrix,
        outputs_matrix,
        algorithm;
        ConvergenceMetric=ConvergenceMetric,
        ConvergenceMetricThreshold=ConvergenceMetricThreshold,
        MaxIter=MaxIter,
        Dampening=Dampening,
        Dampening_With_Input=Dampening_With_Input,
        PrintReports=PrintReports,
        ReportingSigFig=ReportingSigFig,
        ReplaceInvalids=ReplaceInvalids,
        quiet_errors=quiet_errors,
    )
end

# Core implementation: matrix-based (advanced usage)
function fixed_point(
    func::Function,
    inputs::AbstractMatrix{T},
    outputs::AbstractMatrix{<:Number},
    algorithm::FixedPointAlgorithm;
    ConvergenceMetric::Function = (input, output) -> maximum(abs.(output .- input)),
    ConvergenceMetricThreshold::Real = 1e-10,
    MaxIter::Integer = 1000,
    Dampening::Real = 1.0,
    Dampening_With_Input::Bool = false,
    PrintReports::Bool = false,
    ReportingSigFig::Integer = 10,
    ReplaceInvalids::Symbol = :NoAction,
    quiet_errors::Bool = false,
    other_outputs::Union{Missing,NamedTuple} = missing,
) where {T<:Number}
    # Core fixed point iteration algorithm
    # Copy inputs to avoid mutating the original
    Inputs = copy(inputs)
    Outputs = copy(outputs)
    SimpleStartIndex = Integer(size(Outputs)[2])
    if isempty(Outputs)
        if size(Inputs)[2] > 1
            @warn(
                "If you do not give outputs to the function then you can only give one vector of inputs (in a 2d array) to the fixed_pointFunction. So for a function that takes an N dimensional array you should input a Array{Float64}(N,1) array.  As you have input an array of size Array{Float64}(N,k) with k > 1 we have discarded everything but the last column to turn it into a Array{Float64}(N,1) array.\n"
            )
            Inputs = Inputs[:, size(Inputs)[2]]
            Inputs = reshape(Inputs, length(Inputs), 1)
        end
    else
        if size(Inputs) != size(Outputs)
            @warn(
                "If you input a matrix of outputs as well as a matrix of inputs then inputs and outputs must be the same shape. As they differ in this case the last column of the inputs matrix has been taken as the starting point and everything else discarded."
            )
            Inputs = Inputs[:, size(Inputs)[2]]
            Inputs = reshape(Inputs, length(Inputs), 1)
            Outputs = Array{T,2}(undef, size(Inputs)[1], 0)
            SimpleStartIndex = Integer(1)
        end
    end
    LengthOfArray = size(Inputs)[1]
    output_type = promote_type(T, eltype(Outputs))
    final_other_output = other_outputs
    # Do an initial run if no runs have been done:
    if isempty(Outputs)
        ExecutedFunction = execute_function_safely(
            func, Inputs[:, 1]; quiet_errors=quiet_errors
        )
        final_other_output = ExecutedFunction.Other_Output_
        if ExecutedFunction.Error_ != :NoError
            return FixedPointResults(
                Inputs,
                Outputs,
                :InvalidInputOrOutputOfIteration;
                FailedEvaluation_=ExecutedFunction,
                Other_Output=final_other_output,
            )
        end
        output_type = promote_type(typeof(Inputs[1]), typeof(ExecutedFunction.Output_[1]))
        converted_outputs = convert.(Ref(output_type), ExecutedFunction.Output_)
        Outputs = hcat(Outputs, converted_outputs)
        Inputs = convert.(Ref(output_type), Inputs)
    else
        # This ensures that MaxIter refers to max iter excluding any previous passed in results
        MaxIter = MaxIter + size(Outputs)[2]
        # This is to take into account previously passed in simple iterations (without jumps).
        SimpleStartIndex =
            SimpleStartIndex - (size(put_together_without_jumps(Inputs, Outputs))[2])
    end
    # First running through the last column of Inputs to test if we already have a fixed point.
    iter = Integer(size(Outputs)[2])
    # Convergence metrics should always return real numbers, even for complex inputs
    convergence_type = output_type <: Complex ? real(output_type) : output_type
    ConvergenceVector = Array{convergence_type,1}(undef, iter)
    for i in 1:iter
        ConvergenceVector[i] = ConvergenceMetric(Inputs[:, i], Outputs[:, i])
    end
    if ConvergenceVector[iter] < ConvergenceMetricThreshold
        if (PrintReports)
            println(
                "The last column of Inputs matrix is already a fixed point under input convergence metric and convergence threshold",
            )
        end
        return FixedPointResults(
            Inputs,
            Outputs,
            :ReachedConvergenceThreshold;
            ConvergenceVector_=vec(ConvergenceVector),
            Other_Output=final_other_output,
        )
    end
    # Printing a report for initial convergence
    Convergence = ConvergenceVector[iter]
    if (PrintReports)
        println(
            "                                          Algorithm: ",
            lpad(algorithm_name(algorithm), 8),
            ". Iteration: ",
            lpad(iter, 5),
            ". Convergence: ",
            lpad(round(Convergence; sigdigits=ReportingSigFig), ReportingSigFig+4),
            ". Time: ",
            now(),
        )
    end
    iter = iter + 1
    while (Convergence > ConvergenceMetricThreshold) & (iter <= MaxIter)
        # Generate new input using algorithm dispatch
        options = Dict(
            :dampening => Dampening,
            :dampening_with_input => Dampening_With_Input,
            :print_reports => PrintReports,
            :replace_invalids => ReplaceInvalids,
            :simple_start_index => SimpleStartIndex
        )

        NewInputFunctionReturn = fixed_point_new_input(Inputs, Outputs, algorithm, options)

        if PrintReports & !isa(algorithm, Anderson)
            print(lpad("", 42))
        end
        ExecutedFunction = execute_function_safely(
            func, NewInputFunctionReturn; type_check=true, quiet_errors=quiet_errors
        )
        if ExecutedFunction.Error_ != :NoError
            return FixedPointResults(
                Inputs,
                Outputs,
                :InvalidInputOrOutputOfIteration;
                ConvergenceVector_=vec(ConvergenceVector),
                FailedEvaluation_=ExecutedFunction,
                Other_Output=ExecutedFunction.Other_Output_,
            )
        end
        final_other_output = ExecutedFunction.Other_Output_
        Inputs = hcat(Inputs, ExecutedFunction.Input_)
        Outputs = hcat(Outputs, convert(Array{output_type,1}, ExecutedFunction.Output_))
        # Checking and recording convergence
        Convergence = ConvergenceMetric(ExecutedFunction.Input_, ExecutedFunction.Output_)
        ConvergenceVector = vcat(ConvergenceVector, Convergence)
        # Output of report and going to next iteration.
        if (PrintReports)
            println(
                "Algorithm: ",
                lpad(algorithm_name(algorithm), 8),
                ". Iteration: ",
                lpad(iter, 5),
                ". Convergence: ",
                lpad(round(Convergence; sigdigits=ReportingSigFig), ReportingSigFig+4),
                ". Time: ",
                now(),
            )
        end
        iter = iter + 1
    end
    Finish = if (Convergence < ConvergenceMetricThreshold)
        :ReachedConvergenceThreshold
    else
        :ReachedMaxIter
    end
    return FixedPointResults(
        Inputs,
        Outputs,
        Finish;
        ConvergenceVector_=vec(ConvergenceVector),
        Other_Output=final_other_output,
    )
end
