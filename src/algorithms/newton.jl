"""
Newton-based acceleration implementation.
"""

"""
    Newton()

Newton-based acceleration using approximate derivatives.

This method uses finite differences to approximate the derivative and applies
Newton's method to solve f(x) - x = 0.
"""
struct Newton <: FixedPointAlgorithm end

# Algorithm properties
algorithm_name(::Newton) = "Newton"
needs_extrapolation_period(::Newton) = false
get_extrapolation_period(::Newton) = 1
is_polynomial_method(::Newton) = false
is_epsilon_method(::Newton) = false

# Algorithm implementation
"""
Compute the next input using Newton acceleration.
"""
function _compute_proposed_input(inputs, outputs, ::Newton, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]
    simple_start_index = get(options, :simple_start_index, 1)

    if ((completed_iters + simple_start_index) % 2 == 1) & (completed_iters > 1)
        xk1 = inputs[:, (completed_iters - 1)]
        fxk1 = outputs[:, (completed_iters - 1)]
        gxk1 = fxk1 .- xk1
        xk = inputs[:, completed_iters]
        fxk = outputs[:, completed_iters]
        gxk = fxk .- xk
        derivative = (gxk .- gxk1) ./ (xk .- xk1)
        return xk .- (gxk ./ derivative)
    else
        return simple_iterate
    end
end
