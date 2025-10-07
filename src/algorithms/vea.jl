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
algorithm_to_symbol(::VEA) = :VEA

# Algorithm implementation
"""
Compute the next input using Vector Epsilon Algorithm.
"""
function _compute_proposed_input(inputs, outputs, alg::VEA, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    simple_iterates_matrix = put_together_without_jumps(inputs, outputs)
    if (size(simple_iterates_matrix)[2] % alg.extrapolation_period) == 0
        return EpsilonExtrapolation(simple_iterates_matrix, algorithm_to_symbol(alg))
    else
        return simple_iterate
    end
end
