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
    inputs_mat, outputs_mat, output_type, other_output_val, max_iter_adjust, simple_start_index, failure = _initial_evaluation!(
        func, inputs_mat, outputs_mat, options
    )
    !isnothing(failure) && return failure
    max_iterations = options.convergence.max_iterations + max_iter_adjust
    convergence_vector = _compute_initial_convergence(
        inputs_mat, outputs_mat, options, output_type
    )
    if convergence_vector[end] < options.convergence.threshold
        if options.reporting.print_reports
            println(
                "The last column of inputs matrix is already a fixed point under input convergence metric and convergence threshold",
            )
        end
        return FixedPointResults(
            inputs_mat,
            outputs_mat,
            :ReachedConvergenceThreshold;
            convergence_vector=convergence_vector,
            other_output_val=other_output_val,
        )
    end
    _maybe_report_initial(
        options, algorithm, length(convergence_vector), convergence_vector[end]
    )
    result, _, _, _ = _iteration_loop!(
        func,
        inputs_mat,
        outputs_mat,
        algorithm,
        options,
        convergence_vector,
        simple_start_index,
        output_type,
        other_output_val,
        max_iterations,
    )
    return result
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
    # Returns (updated_inputs, updated_outputs, output_type, other_output, max_iterations_adjust, simple_start_index_adjust, failure_result|nothing)
    simple_start_index = size(outputs_mat, 2)
    max_iterations_adjust = 0
    simple_start_index_adjust = 0
    other_output_val = missing
    output_type = promote_type(T, eltype(outputs_mat))
    if isempty(outputs_mat)
        executed_func = execute_function_safely(
            func, inputs_mat[:, 1]; quiet_errors=options.stability.quiet_errors
        )
        other_output_val = executed_func.other_output
        if executed_func.error != :NoError
            return (
                inputs_mat,
                outputs_mat,
                output_type,
                other_output_val,
                max_iterations_adjust,
                simple_start_index_adjust,
                FixedPointResults(
                    inputs_mat,
                    outputs_mat,
                    :InvalidInputOrOutputOfIteration;
                    failed_evaluation=executed_func,
                    other_output_val=other_output_val,
                ),
            )
        end
        output_type = promote_type(typeof(inputs_mat[1]), typeof(executed_func.output[1]))
        converted_outputs = convert.(Ref(output_type), executed_func.output)
        outputs_mat = hcat(outputs_mat, converted_outputs)
        inputs_mat = convert.(Ref(output_type), inputs_mat)
    else
        max_iterations_adjust = size(outputs_mat, 2)
        simple_start_index_adjust =
            -size(put_together_without_jumps(inputs_mat, outputs_mat), 2)
    end
    return (
        inputs_mat,
        outputs_mat,
        output_type,
        other_output_val,
        max_iterations_adjust,
        simple_start_index + simple_start_index_adjust,
        nothing,
    )
end

function _compute_initial_convergence(
    inputs_mat, outputs_mat, options::FixedPointOptions, output_type
)
    iteration = size(outputs_mat, 2)
    convergence_type = output_type <: Complex ? real(output_type) : output_type
    convergence_vector = Vector{convergence_type}(undef, iteration)
    metric = options.convergence.metric
    for i in 1:iteration
        convergence_vector[i] = metric(inputs_mat[:, i], outputs_mat[:, i])
    end
    return convergence_vector
end

function _maybe_report_initial(
    options::FixedPointOptions, algorithm::FixedPointAlgorithm, iteration::Int, convergence
)
    if options.reporting.print_reports
        println(
            "                                          Algorithm: ",
            lpad(algorithm_name(algorithm), 8),
            ". Iteration: ",
            lpad(iteration, 5),
            ". Convergence: ",
            lpad(
                round(convergence; sigdigits=options.reporting.reporting_sig_figs),
                options.reporting.reporting_sig_figs + 4,
            ),
            ". Time: ",
            now(),
        )
    end
end

function _iteration_loop!(
    func::Function,
    inputs_mat::Matrix{T},
    outputs_mat::Matrix{T},
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions,
    convergence_vector::Vector{<:Real},
    simple_start_index::Int,
    output_type,
    other_output_val,
    max_iterations::Int,
) where {T<:Number}
    convergence = convergence_vector[end]
    iteration = length(convergence_vector) + 1
    while convergence > options.convergence.threshold && iteration <= max_iterations
        new_input = fixed_point_new_input(
            inputs_mat, outputs_mat, algorithm, options, simple_start_index
        )
        if options.reporting.print_reports && !isa(algorithm, Anderson)
            print(lpad("", 42))
        end
        executed_func = execute_function_safely(
            func, new_input; type_check=true, quiet_errors=options.stability.quiet_errors
        )
        if executed_func.error != :NoError
            return (
                FixedPointResults(
                    inputs_mat,
                    outputs_mat,
                    :InvalidInputOrOutputOfIteration;
                    convergence_vector=convergence_vector,
                    failed_evaluation=executed_func,
                    other_output_val=executed_func.other_output,
                ),
                inputs_mat,
                outputs_mat,
                other_output_val,
            )
        end
        other_output_val = executed_func.other_output
        inputs_mat = hcat(inputs_mat, executed_func.input)
        outputs_mat = hcat(outputs_mat, convert(Vector{output_type}, executed_func.output))
        convergence = options.convergence.metric(executed_func.input, executed_func.output)
        push!(convergence_vector, convergence)
        if options.reporting.print_reports
            println(
                "Algorithm: ",
                lpad(algorithm_name(algorithm), 8),
                ". Iteration: ",
                lpad(iteration, 5),
                ". Convergence: ",
                lpad(
                    round(convergence; sigdigits=options.reporting.reporting_sig_figs),
                    options.reporting.reporting_sig_figs + 4,
                ),
                ". Time: ",
                now(),
            )
        end
        iteration += 1
    end
    final_status = if convergence < options.convergence.threshold
        :ReachedConvergenceThreshold
    else
        :ReachedMaxIter
    end
    return (
        FixedPointResults(
            inputs_mat,
            outputs_mat,
            final_status;
            convergence_vector=convergence_vector,
            other_output_val=other_output_val,
        ),
        inputs_mat,
        outputs_mat,
        other_output_val,
    )
end
