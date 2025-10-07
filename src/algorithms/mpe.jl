"""
Minimal Polynomial Extrapolation implementation.
"""

"""
    MPE(; extrapolation_period=7)

Minimal Polynomial Extrapolation.

MPE takes a linear combination of previous sequential iterates to accelerate convergence.
Extrapolation is performed every `extrapolation_period` iterations.

# Parameters
- `extrapolation_period::Int`: Number of simple iterations before extrapolation (default: 7)

# References
Cabay, S., and L.W. Jackson (1976). "A Polynomial Extrapolation Method for Finding
Limits and Antilimits of Vector Sequences." SIAM Journal of Numerical Analysis 13(5): 734-52.
"""
struct MPE <: FixedPointAlgorithm
    extrapolation_period::Int

    function MPE(; extrapolation_period::Int=7)
        extrapolation_period >= 1 ||
            throw(ArgumentError("extrapolation_period must be at least 1"))
        new(extrapolation_period)
    end
end

# Algorithm properties
algorithm_name(::MPE) = "MPE"
needs_extrapolation_period(::MPE) = true
get_extrapolation_period(alg::MPE) = alg.extrapolation_period
is_polynomial_method(::MPE) = true
is_epsilon_method(::MPE) = false
algorithm_to_symbol(::MPE) = :MPE

# Algorithm implementation
"""
Compute the next input using Minimal Polynomial Extrapolation.
"""
function _compute_proposed_input(inputs, outputs, alg::MPE, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    simple_iterates_matrix = put_together_without_jumps(inputs, outputs)
    if (size(simple_iterates_matrix)[2] % alg.extrapolation_period == 0)
        return PolynomialExtrapolation(simple_iterates_matrix, algorithm_to_symbol(alg))
    else
        return simple_iterate
    end
end
