"""
    fixed_point(func, initial_guess, algorithm, options::FixedPointOptions = default_options())

Find the fixed point of a function using various acceleration algorithms.

# Arguments
- `func::Function`: The function for which a fixed point is sought
- `initial_guess`: Initial guess for the fixed point (scalar, vector, or matrix)
- `algorithm::FixedPointAlgorithm`: Algorithm instance (e.g., Anderson(), Simple())
- `options::FixedPointOptions`: Configuration options (optional, defaults to default_options())

# Examples
```julia
# Use default options
result = fixed_point(f, 0.3, Aitken())

# Use preset options
result = fixed_point(f, 0.3, Aitken(), robust_options())

# Custom configuration
opts = FixedPointOptions(threshold=1e-12, print_reports=true)
result = fixed_point(g, [0.3, 900.0], Anderson(), opts)
```
"""
function fixed_point(
    func::Function,
    initial_guess::Number,
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions = default_options()
)
    inputs_matrix = Array{typeof(initial_guess),2}(undef, 1, 1)
    inputs_matrix[1, 1] = initial_guess
    outputs_matrix = Array{typeof(initial_guess),2}(undef, 1, 0)
    return fixed_point(
        func,
        inputs_matrix,
        outputs_matrix,
        algorithm,
        options
    )
end

function fixed_point(
    func::Function,
    initial_guess::AbstractVector{T},
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions = default_options()
) where {T<:Number}
    inputs_matrix = reshape(initial_guess, length(initial_guess), 1)
    outputs_matrix = Array{T,2}(undef, size(inputs_matrix)[1], 0)
    return fixed_point(
        func,
        inputs_matrix,
        outputs_matrix,
        algorithm,
        options
    )
end

"""
    fixed_point(func, previous_results, algorithm, options::FixedPointOptions = default_options())

Continue from previous results.
"""
function fixed_point(
    func::Function,
    previous_results::FixedPointResults,
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions = default_options()
)
    return fixed_point(
        func,
        previous_results.Inputs_,
        previous_results.Outputs_,
        algorithm,
        options
    )
end

"""
    fixed_point(func, inputs_matrix, outputs_matrix, algorithm, options::FixedPointOptions = default_options())

Core matrix interface for advanced usage.
"""
function fixed_point(
    func::Function,
    inputs_matrix::AbstractMatrix,
    outputs_matrix::AbstractMatrix,
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions = default_options()
)
    return _fixed_point_core(
        func,
        inputs_matrix,
        outputs_matrix,
        algorithm;
        ConvergenceMetric = options.convergence.metric,
        ConvergenceMetricThreshold = options.convergence.threshold,
        MaxIter = options.convergence.max_iterations,
        Dampening = options.stability.dampening,
        Dampening_With_Input = options.stability.dampening_with_input,
        ReplaceInvalids = options.stability.replace_invalids,
        quiet_errors = options.stability.quiet_errors,
        PrintReports = options.reporting.print_reports,
        ReportingSigFig = options.reporting.reporting_sig_figs
    )
end

# Core implementation: matrix-based (advanced usage)
function _fixed_point_core(
    func::Function,
    inputs::AbstractMatrix{T},
    outputs::AbstractMatrix{<:Number},
    algorithm::FixedPointAlgorithm;
    ConvergenceMetric::Function=(input, output) -> maximum(abs.(output .- input)),
    ConvergenceMetricThreshold::Real=1e-10,
    MaxIter::Integer=1000,
    Dampening::Real=1.0,
    Dampening_With_Input::Bool=false,
    PrintReports::Bool=false,
    ReportingSigFig::Integer=10,
    ReplaceInvalids::Symbol=:NoAction,
    quiet_errors::Bool=false,
) where {T<:Number}
    # Core fixed point iteration algorithm
    # Copy inputs to avoid mutating the original
    Inputs = copy(inputs)
    Outputs = copy(outputs)
    SimpleStartIndex = Integer(size(Outputs)[2])
    other_outputs = missing  # For compatibility
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
            :simple_start_index => SimpleStartIndex,
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
