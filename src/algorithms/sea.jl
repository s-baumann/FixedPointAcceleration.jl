"""
Scalar Epsilon Algorithm implementation.
"""

"""
    SEA(; extrapolation_period=7)

Scalar Epsilon Algorithm.

SEA applies the epsilon algorithm elementwise to vector problems.
Extrapolation is performed every `extrapolation_period` iterations.

# Parameters
- `extrapolation_period::Int`: Number of simple iterations before extrapolation (default: 7)
                              Should be even for optimal performance.

# References
Wynn, P. (1956). "On a Device for Computing the e_m(S_n) Transformation."
Mathematical Tables and Other Aids to Computation 10(54): 91-96.
"""
struct SEA <: FixedPointAlgorithm
    extrapolation_period::Int

    function SEA(; extrapolation_period::Int=7)
        extrapolation_period >= 1 ||
            throw(ArgumentError("extrapolation_period must be at least 1"))
        new(extrapolation_period)
    end
end

# Algorithm properties
algorithm_name(::SEA) = "SEA"
needs_extrapolation_period(::SEA) = true
get_extrapolation_period(alg::SEA) = alg.extrapolation_period
is_polynomial_method(::SEA) = false
is_epsilon_method(::SEA) = true


# Algorithm implementation
"""
Compute the next input using Scalar Epsilon Algorithm.
"""
function _compute_proposed_input(inputs, outputs, alg::SEA, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    simple_iterates_matrix = put_together_without_jumps(inputs, outputs)
    if (size(simple_iterates_matrix)[2] % alg.extrapolation_period) == 0
        return _sea_epsilon_extrapolation(simple_iterates_matrix)
    else
        return simple_iterate
    end
end

"""
Perform Scalar Epsilon Algorithm extrapolation on a matrix of iterates.
"""
function _sea_epsilon_extrapolation(iterates::AbstractArray{R,2}) where {R<:Number}
    # Handle single column case
    if size(iterates)[2] == 1
        return iterates[:, 1]
    end

    # Ensure odd number of columns for epsilon algorithm
    if size(iterates)[2] % 2 == 0
        iterates = iterates[:, 2:end]
    end

    mat = iterates
    rows_of_matrix = size(mat)[1]
    total_columns = size(mat)[2]
    previous_matrix = zeros(rows_of_matrix, total_columns - 1)

    for matrix_column in reverse(2:total_columns)
        diff_matrix = mat[:, 2:matrix_column] .- mat[:, 1:(matrix_column - 1)]
        new_matrix = previous_matrix + (1 ./ diff_matrix)  # SEA uses elementwise inverse
        previous_matrix = mat[:, 2:(matrix_column - 1)]
        mat = new_matrix
    end

    # Handle NaN/missing values by trying with fewer columns
    if eltype(mat) <: Complex
        if any(isnan.(real.(mat))) || any(isnan.(imag.(mat))) || any(ismissing.(mat))
            return _sea_epsilon_extrapolation(iterates[:, 3:end])
        end
    else
        if any(isnan.(mat)) || any(ismissing.(mat))
            return _sea_epsilon_extrapolation(iterates[:, 3:end])
        end
    end

    return mat[:, 1]
end
