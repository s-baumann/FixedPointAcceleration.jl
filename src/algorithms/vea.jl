"""
Vector Epsilon Algorithm implementation.
"""

"""
    VEA(; extrapolation_period=7)

Vector Epsilon Algorithm.

VEA uses the epsilon algorithm with Moore-Penrose pseudoinverse for vector problems.
Extrapolation is performed every `extrapolation_period` iterations.

# Parameters
- `extrapolation_period::Int`: Number of simple iterations before extrapolation (default: 7)
                              Should be even for optimal performance.

# References
Wynn, P. (1956). "On a Device for Computing the e_m(S_n) Transformation."
Mathematical Tables and Other Aids to Computation 10(54): 91-96.
"""
struct VEA <: FixedPointAlgorithm
    extrapolation_period::Int

    function VEA(; extrapolation_period::Int=7)
        extrapolation_period >= 1 ||
            throw(ArgumentError("extrapolation_period must be at least 1"))
        new(extrapolation_period)
    end
end

# Algorithm properties
algorithm_name(::VEA) = "VEA"
needs_extrapolation_period(::VEA) = true
get_extrapolation_period(alg::VEA) = alg.extrapolation_period
is_polynomial_method(::VEA) = false
is_epsilon_method(::VEA) = true


# Algorithm implementation
"""
Compute the next input using Vector Epsilon Algorithm.
"""
function _compute_proposed_input(inputs, outputs, alg::VEA, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    simple_iterates_matrix = put_together_without_jumps(inputs, outputs)
    if (size(simple_iterates_matrix)[2] % alg.extrapolation_period) == 0
        return _vea_epsilon_extrapolation(simple_iterates_matrix)
    else
        return simple_iterate
    end
end

"""
Perform Vector Epsilon Algorithm extrapolation on a matrix of iterates.
"""
function _vea_epsilon_extrapolation(iterates::AbstractArray{R,2}) where {R<:Number}
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
        new_matrix = previous_matrix + _vea_vector_of_inverses(diff_matrix)
        previous_matrix = mat[:, 2:(matrix_column - 1)]
        mat = new_matrix
    end

    # Handle NaN/missing values by trying with fewer columns
    if eltype(mat) <: Complex
        if any(isnan.(real.(mat))) || any(isnan.(imag.(mat))) || any(ismissing.(mat))
            return _vea_epsilon_extrapolation(iterates[:, 3:end])
        end
    else
        if any(isnan.(mat)) || any(ismissing.(mat))
            return _vea_epsilon_extrapolation(iterates[:, 3:end])
        end
    end

    return mat[:, 1]
end

"""
Compute vector of Moore-Penrose pseudoinverses for VEA.
"""
function _vea_vector_of_inverses(difference_matrix::AbstractArray{T,2}) where {T<:Number}
    if size(difference_matrix)[1] == 1
        return 1 ./ difference_matrix
    else
        invs = transpose(pinv(difference_matrix[:, 1]))
        if size(difference_matrix)[2] < 2
            return invs
        end
        for i in 2:size(difference_matrix)[2]
            invs = hcat(invs, transpose(pinv(difference_matrix[:, i])))
        end
        return invs
    end
end
