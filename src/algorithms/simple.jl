"""
Simple fixed point iteration (no acceleration).
"""

"""
    Simple()

Simple fixed point iteration: x_{n+1} = f(x_n)

This is the most basic fixed point iteration method with no acceleration.
It's guaranteed to converge for contractive functions but may be slow.
"""
struct Simple <: FixedPointAlgorithm end

# Algorithm properties
algorithm_name(::Simple) = "Simple"
needs_extrapolation_period(::Simple) = false
get_extrapolation_period(::Simple) = 1
is_polynomial_method(::Simple) = false
is_epsilon_method(::Simple) = false

# Algorithm implementation
"""
Compute the next input using simple iteration.
"""
function _compute_proposed_input(
    inputs, outputs, ::Simple, options::FixedPointOptions, simple_start_index::Int
)
    return outputs[:, end]
end
