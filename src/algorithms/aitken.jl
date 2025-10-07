"""
Aitken's Δ² acceleration method implementation.
"""

"""
    Aitken()

Aitken's Δ² acceleration method.

This scalar acceleration method (applied elementwise to vectors) assumes the sequence
converges at a constant rate and uses three consecutive iterates to estimate the limit.

# References
Aitken, A.C. (1926). "On Bernoulli's Numerical Solution of Algebraic Equations."
Proceedings of the Royal Society of Edinburgh 46: 289-305.
"""
struct Aitken <: FixedPointAlgorithm end

# Algorithm properties
algorithm_name(::Aitken) = "Aitken"
needs_extrapolation_period(::Aitken) = false
get_extrapolation_period(::Aitken) = 1
is_polynomial_method(::Aitken) = false
is_epsilon_method(::Aitken) = false
algorithm_to_symbol(::Aitken) = :Aitken

# Algorithm implementation
"""
Compute the next input using Aitken acceleration.
"""
function _compute_proposed_input(inputs, outputs, ::Aitken, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]
    simple_start_index = get(options, :simple_start_index, 1)

    if ((completed_iters + simple_start_index) % 3) == 0
        x = inputs[:, (completed_iters - 1)]
        fx = outputs[:, (completed_iters - 1)]
        ffx = outputs[:, completed_iters]
        return x .- ((fx .- x) .^ 2 ./ (ffx .- 2 .* fx .+ x))
    else
        return simple_iterate
    end
end
