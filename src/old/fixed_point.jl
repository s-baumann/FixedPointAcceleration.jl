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
    inputs::AbstractMatrix{T},
    outputs::AbstractMatrix{<:Number},
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions,
) where {T<:Number}
    inputs_mat, outputs_mat, simple_start_index = _prepare_inputs_outputs(inputs, outputs)
    state = _initial_evaluation!(func, inputs_mat, outputs_mat, options)

    # Check if there was a failure during initial evaluation
    if state.other_output_val isa FixedPointResults
        return state.other_output_val
    end

    # Check if already converged
    if state.convergence_vector[end] < options.threshold
        if options.print_reports
            @warn "The last column of inputs matrix is already a fixed point under input convergence metric and convergence threshold"
        end
        return FixedPointResults(state, :ReachedConvergenceThreshold)
    end

    _maybe_report_initial(
        options, algorithm, length(state.convergence_vector), state.convergence_vector[end]
    )

    return _iteration_loop!(func, algorithm, options, state)
end

# === Refactored internal pipeline ===

function _prepare_inputs_outputs(
    inputs::AbstractMatrix{T}, outputs::AbstractMatrix{<:Number}
) where {T<:Number}
    inputs_mat = copy(inputs)
    outputs_mat = copy(outputs)
    simple_start_index = size(outputs_mat, 2)
    if isempty(outputs_mat)
        if size(inputs_mat, 2) > 1
            @warn "If you do not give outputs to the function then you can only give one vector of inputs (in a 2d array) to the fixed_pointFunction. So for a function that takes an N dimensional array you should input a Array{Float64}(N,1) array.  As you have input an array of size Array{Float64}(N,k) with k > 1 we have discarded everything but the last column to turn it into a Array{Float64}(N,1) array.\n"
            inputs_mat = inputs_mat[:, end:end]
        end
    else
        if size(inputs_mat) != size(outputs_mat)
            @warn "If you input a matrix of outputs as well as a matrix of inputs then inputs and outputs must be the same shape. As they differ in this case the last column of the inputs matrix has been taken as the starting point and everything else discarded."
            inputs_mat = inputs_mat[:, end:end]
            outputs_mat = Matrix{T}(undef, size(inputs_mat, 1), 0)
            simple_start_index = 1
        end
    end
    return inputs_mat, outputs_mat, simple_start_index
end

function _initial_evaluation!(
    func::Function,
    inputs_mat::Matrix{T},
    outputs_mat::Matrix{R},
    options::FixedPointOptions,
) where {T<:Number,R<:Number}
    # Always returns _IterationState - type stable
    simple_start_index = size(outputs_mat, 2)
    max_iterations_adjust = 0
    simple_start_index_adjust = 0
    other_output_val = missing
    output_type = promote_type(T, eltype(outputs_mat))
    failure_result = nothing

    if isempty(outputs_mat)
        executed_func = execute_function_safely(
            func, inputs_mat[:, 1]; quiet_errors=options.quiet_errors
        )
        other_output_val = executed_func.other_output
        if executed_func.error != :NoError
            # Store failure result to be handled by caller
            failure_result = FixedPointResults(
                inputs_mat,
                outputs_mat,
                :InvalidInputOrOutputOfIteration;
                failed_evaluation=executed_func,
                other_output_val=other_output_val,
            )
            # Create a dummy convergence vector to avoid errors
            convergence_vector = [Inf]
        else
            output_type = promote_type(
                typeof(inputs_mat[1]), typeof(executed_func.output[1])
            )
            converted_outputs = convert.(Ref(output_type), executed_func.output)
            outputs_mat = hcat(outputs_mat, converted_outputs)
            inputs_mat = convert.(Ref(output_type), inputs_mat)
            # Compute initial convergence
            convergence_vector = _compute_initial_convergence(
                inputs_mat, outputs_mat, options, output_type
            )
        end
    else
        max_iterations_adjust = size(outputs_mat, 2)
        simple_start_index_adjust =
            -size(put_together_without_jumps(inputs_mat, outputs_mat), 2)
        # Compute initial convergence
        convergence_vector = _compute_initial_convergence(
            inputs_mat, outputs_mat, options, output_type
        )
    end

    max_iterations = options.max_iterations + max_iterations_adjust

    # Construct and return the iteration state
    state = _IterationState(
        inputs_mat,
        outputs_mat,
        convergence_vector,
        simple_start_index + simple_start_index_adjust,
        other_output_val,
        max_iterations,
    )

    # Store failure result in the state for the caller to check
    if !isnothing(failure_result)
        state.other_output_val = failure_result
    end

    return state
end

function _compute_initial_convergence(
    inputs_mat, outputs_mat, options::FixedPointOptions, output_type
)
    iteration = size(outputs_mat, 2)
    convergence_type = output_type <: Complex ? real(output_type) : output_type
    convergence_vector = Vector{convergence_type}(undef, iteration)
    metric = options.metric
    for i in 1:iteration
        convergence_vector[i] = metric(inputs_mat[:, i], outputs_mat[:, i])
    end
    return convergence_vector
end

function _maybe_report_initial(
    options::FixedPointOptions, algorithm::FixedPointAlgorithm, iteration::Int, convergence
)
    if options.print_reports
        println(
            "                                          Algorithm: ",
            lpad(algorithm_name(algorithm), 8),
            ". Iteration: ",
            lpad(iteration, 5),
            ". Convergence: ",
            lpad(
                round(convergence; sigdigits=options.reporting_sig_figs),
                options.reporting_sig_figs + 4,
            ),
            ". Time: ",
            now(),
        )
    end
end

function _iteration_loop!(
    func::Function,
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions,
    state::_IterationState,
)
    convergence = state.convergence_vector[end]
    iteration = length(state.convergence_vector) + 1
    while convergence > options.threshold && iteration <= state.max_iterations
        new_input = fixed_point_new_input(
            state.inputs, state.outputs, algorithm, options, state.simple_start_index
        )
        if options.print_reports && !isa(algorithm, Anderson)
            print(lpad("", 42))
        end
        executed_func = execute_function_safely(
            func, new_input; type_check=true, quiet_errors=options.quiet_errors
        )
        if executed_func.error != :NoError
            return FixedPointResults(
                state.inputs,
                state.outputs,
                :InvalidInputOrOutputOfIteration;
                convergence_vector=state.convergence_vector,
                failed_evaluation=executed_func,
                other_output_val=executed_func.other_output,
            )
        end
        state.other_output_val = executed_func.other_output
        state.inputs = hcat(state.inputs, executed_func.input)
        col_type = eltype(state.inputs)
        state.outputs = hcat(state.outputs, convert(Vector{col_type}, executed_func.output))
        convergence = options.metric(executed_func.input, executed_func.output)
        push!(state.convergence_vector, convergence)

        _maybe_report_initial(options, algorithm, iteration, convergence)
        iteration += 1
    end
    final_status =
        convergence < options.threshold ? :ReachedConvergenceThreshold : :ReachedMaxIter
    return FixedPointResults(
        state.inputs,
        state.outputs,
        final_status;
        convergence_vector=state.convergence_vector,
        other_output_val=state.other_output_val,
    )
end
