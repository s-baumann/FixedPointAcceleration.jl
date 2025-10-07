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
    options::FixedPointOptions=default_options(),
)
    T = typeof(initial_guess)
    inputs_matrix = Matrix{T}(undef, 1, 1)
    inputs_matrix[1, 1] = initial_guess
    outputs_matrix = Matrix{T}(undef, 1, 0)
    return fixed_point(func, inputs_matrix, outputs_matrix, algorithm, options)
end

function fixed_point(
    func::Function,
    initial_guess::AbstractVector{T},
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions=default_options(),
) where {T<:Number}
    inputs_matrix = reshape(initial_guess, length(initial_guess), 1)
    outputs_matrix = Matrix{T}(undef, size(inputs_matrix, 1), 0)
    return fixed_point(func, inputs_matrix, outputs_matrix, algorithm, options)
end

"""
    fixed_point(func, previous_results, algorithm, options::FixedPointOptions = default_options())

Continue from previous results.
"""
function fixed_point(
    func::Function,
    previous_results::FixedPointResults,
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions=default_options(),
)
    return fixed_point(
        func, previous_results.inputs, previous_results.outputs, algorithm, options
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
    options::FixedPointOptions=default_options(),
)
    return _fixed_point_core(
        func,
        inputs_matrix,
        outputs_matrix,
        algorithm;
        convergence_metric=options.convergence.metric,
        convergence_threshold=options.convergence.threshold,
        max_iterations=options.convergence.max_iterations,
        dampening=options.stability.dampening,
        dampening_with_input=options.stability.dampening_with_input,
        replace_invalids=options.stability.replace_invalids,
        quiet_errors=options.stability.quiet_errors,
        print_reports=options.reporting.print_reports,
        reporting_sig_figs=options.reporting.reporting_sig_figs,
    )
end

# Core implementation: matrix-based (advanced usage)
function _fixed_point_core(
    func::Function,
    inputs::AbstractMatrix{T},
    outputs::AbstractMatrix{<:Number},
    algorithm::FixedPointAlgorithm;
    convergence_metric::Function=(input, output) -> maximum(abs.(output .- input)),
    convergence_threshold::Real=1e-10,
    max_iterations::Integer=1000,
    dampening::Real=1.0,
    dampening_with_input::Bool=false,
    print_reports::Bool=false,
    reporting_sig_figs::Integer=10,
    replace_invalids::Symbol=:NoAction,
    quiet_errors::Bool=false,
) where {T<:Number}
    # Core fixed point iteration algorithm
    # Copy inputs to avoid mutating the original
    inputs_mat = copy(inputs)
    outputs_mat = copy(outputs)
    simple_start_index = size(outputs_mat, 2)
    other_outputs = missing  # For compatibility

    if isempty(outputs_mat)
        if size(inputs_mat, 2) > 1
            @warn(
                "If you do not give outputs to the function then you can only give one vector of inputs (in a 2d array) to the fixed_pointFunction. So for a function that takes an N dimensional array you should input a Array{Float64}(N,1) array.  As you have input an array of size Array{Float64}(N,k) with k > 1 we have discarded everything but the last column to turn it into a Array{Float64}(N,1) array.\n"
            )
            inputs_mat = inputs_mat[:, end:end]
        end
    else
        if size(inputs_mat) != size(outputs_mat)
            @warn(
                "If you input a matrix of outputs as well as a matrix of inputs then inputs and outputs must be the same shape. As they differ in this case the last column of the inputs matrix has been taken as the starting point and everything else discarded."
            )
            inputs_mat = inputs_mat[:, end:end]
            outputs_mat = Matrix{T}(undef, size(inputs_mat, 1), 0)
            simple_start_index = 1
        end
    end

    array_length = size(inputs_mat, 1)
    output_type = promote_type(T, eltype(outputs_mat))
    final_other_output = other_outputs
    # Do an initial run if no runs have been done
    if isempty(outputs_mat)
        executed_func = execute_function_safely(
            func, inputs_mat[:, 1]; quiet_errors=quiet_errors
        )
        final_other_output = executed_func.other_output
        if executed_func.error != :NoError
            return FixedPointResults(
                inputs_mat,
                outputs_mat,
                :InvalidInputOrOutputOfIteration;
                failed_evaluation=executed_func,
                other_output_val=final_other_output,
            )
        end
        output_type = promote_type(typeof(inputs_mat[1]), typeof(executed_func.output[1]))
        converted_outputs = convert.(Ref(output_type), executed_func.output)
        outputs_mat = hcat(outputs_mat, converted_outputs)
        inputs_mat = convert.(Ref(output_type), inputs_mat)
    else
        # This ensures that max_iterations refers to max iter excluding any previous passed in results
        max_iterations = max_iterations + size(outputs_mat, 2)
        # This is to take into account previously passed in simple iterations (without jumps).
        simple_start_index =
            simple_start_index -
            size(put_together_without_jumps(inputs_mat, outputs_mat), 2)
    end
    # First running through the last column of inputs to test if we already have a fixed point
    iteration = size(outputs_mat, 2)
    # Convergence metrics should always return real numbers, even for complex inputs
    convergence_type = output_type <: Complex ? real(output_type) : output_type
    convergence_vector = Vector{convergence_type}(undef, iteration)
    for i in 1:iteration
        convergence_vector[i] = convergence_metric(inputs_mat[:, i], outputs_mat[:, i])
    end

    if convergence_vector[iteration] < convergence_threshold
        if print_reports
            println(
                "The last column of inputs matrix is already a fixed point under input convergence metric and convergence threshold",
            )
        end
        return FixedPointResults(
            inputs_mat,
            outputs_mat,
            :ReachedConvergenceThreshold;
            convergence_vector=convergence_vector,
            other_output_val=final_other_output,
        )
    end
    # Printing a report for initial convergence
    convergence = convergence_vector[iteration]
    if print_reports
        println(
            "                                          Algorithm: ",
            lpad(algorithm_name(algorithm), 8),
            ". Iteration: ",
            lpad(iteration, 5),
            ". Convergence: ",
            lpad(round(convergence; sigdigits=reporting_sig_figs), reporting_sig_figs+4),
            ". Time: ",
            now(),
        )
    end
    iteration += 1

    while convergence > convergence_threshold && iteration <= max_iterations
        # Generate new input using algorithm dispatch
        algorithm_options = Dict(
            :dampening => dampening,
            :dampening_with_input => dampening_with_input,
            :print_reports => print_reports,
            :replace_invalids => replace_invalids,
            :simple_start_index => simple_start_index,
        )

        new_input = fixed_point_new_input(
            inputs_mat, outputs_mat, algorithm, algorithm_options
        )

        if print_reports && !isa(algorithm, Anderson)
            print(lpad("", 42))
        end

        executed_func = execute_function_safely(
            func, new_input; type_check=true, quiet_errors=quiet_errors
        )

        if executed_func.error != :NoError
            return FixedPointResults(
                inputs_mat,
                outputs_mat,
                :InvalidInputOrOutputOfIteration;
                convergence_vector=convergence_vector,
                failed_evaluation=executed_func,
                other_output_val=executed_func.other_output,
            )
        end

        final_other_output = executed_func.other_output
        inputs_mat = hcat(inputs_mat, executed_func.input)
        outputs_mat = hcat(outputs_mat, convert(Vector{output_type}, executed_func.output))

        # Check and record convergence
        convergence = convergence_metric(executed_func.input, executed_func.output)
        push!(convergence_vector, convergence)

        # Output report and go to next iteration
        if print_reports
            println(
                "Algorithm: ",
                lpad(algorithm_name(algorithm), 8),
                ". Iteration: ",
                lpad(iteration, 5),
                ". Convergence: ",
                lpad(
                    round(convergence; sigdigits=reporting_sig_figs), reporting_sig_figs+4
                ),
                ". Time: ",
                now(),
            )
        end
        iteration += 1
    end

    # Determine final status
    final_status = if convergence < convergence_threshold
        :ReachedConvergenceThreshold
    else
        :ReachedMaxIter
    end

    return FixedPointResults(
        inputs_mat,
        outputs_mat,
        final_status;
        convergence_vector=convergence_vector,
        other_output_val=final_other_output,
    )
end
