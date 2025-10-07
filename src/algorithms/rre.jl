"""
Reduced Rank Extrapolation implementation.
"""

"""
    RRE(; extrapolation_period=7)

Reduced Rank Extrapolation.

RRE is similar to MPE but uses a different linear combination approach.
Extrapolation is performed every `extrapolation_period` iterations.

# Parameters
- `extrapolation_period::Int`: Number of simple iterations before extrapolation (default: 7)

# References
Eddy, R.P. (1979). "Extrapolating to the Limit of a Vector Sequence."
Information Linkage Between Applied Mathematics and Industry: 387-96.
"""
struct RRE <: FixedPointAlgorithm
    extrapolation_period::Int

    function RRE(; extrapolation_period::Int=7)
        extrapolation_period >= 1 ||
            throw(ArgumentError("extrapolation_period must be at least 1"))
        new(extrapolation_period)
    end
end

# Algorithm properties
algorithm_name(::RRE) = "RRE"
needs_extrapolation_period(::RRE) = true
get_extrapolation_period(alg::RRE) = alg.extrapolation_period
is_polynomial_method(::RRE) = true
is_epsilon_method(::RRE) = false


# Algorithm implementation
"""
Compute the next input using Reduced Rank Extrapolation.
"""
function _compute_proposed_input(inputs, outputs, alg::RRE, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    simple_iterates_matrix = put_together_without_jumps(inputs, outputs)
    if (size(simple_iterates_matrix)[2] % alg.extrapolation_period == 0)
        return _rre_extrapolation(simple_iterates_matrix)
    else
        return simple_iterate
    end
end

"""
Perform Reduced Rank Extrapolation on a matrix of iterates.
"""
function _rre_extrapolation(iterates::AbstractArray{R,2}) where {R<:Number}
    total_columns = size(iterates)[2]
    first_column = iterates[:, 1]
    differences = iterates[:, 2:total_columns] - iterates[:, 1:(total_columns - 1)]
    second_differences = differences[:, 2:(total_columns - 1)] - differences[:, 1:(total_columns - 2)]
    first_difference = differences[:, 1]
    differences = differences[:, 1:(total_columns - 2)]
    inverse_second_differences = pinv(second_differences)
    return first_column - ((differences * inverse_second_differences) * first_difference)
end
