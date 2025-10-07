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
algorithm_to_symbol(::RRE) = :RRE

# Algorithm implementation
"""
Compute the next input using Reduced Rank Extrapolation.
"""
function _compute_proposed_input(inputs, outputs, alg::RRE, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    simple_iterates_matrix = put_together_without_jumps(inputs, outputs)
    if (size(simple_iterates_matrix)[2] % alg.extrapolation_period == 0)
        return PolynomialExtrapolation(simple_iterates_matrix, algorithm_to_symbol(alg))
    else
        return simple_iterate
    end
end
