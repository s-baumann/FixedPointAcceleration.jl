"""
Internal function to compute the raw proposed input (before post-processing).
This is where the algorithm-specific dispatch happens.
The actual implementations are in the individual algorithm files.
"""
function _compute_proposed_input end

"""
    fixed_point_new_input(inputs, outputs, algorithm::FixedPointAlgorithm, options::Dict)

Modern dispatch-based interface for computing the next input in fixed point iteration.

### Inputs
* `inputs` - N×A matrix of previous inputs
* `outputs` - N×A matrix of corresponding outputs
* `algorithm` - Algorithm instance (e.g., Anderson(), Simple(), etc.)
* `options` - Dictionary containing algorithm options

### Returns
* Vector of the next guess for the fixed point
"""
function fixed_point_new_input(
    inputs::AbstractArray{T,2},
    outputs::AbstractArray{T,2},
    algorithm::FixedPointAlgorithm,
    options::Dict,
) where {T<:Number}
    # Get the raw proposed input from the algorithm-specific implementation
    proposed_input = _compute_proposed_input(inputs, outputs, algorithm, options)

    # Apply dampening and replacement strategies
    return _apply_post_processing(inputs, outputs, proposed_input, options)
end

"""
Apply post-processing: replacement strategies and dampening.
"""
function _apply_post_processing(
    inputs::AbstractArray{T,2},
    outputs::AbstractArray{T,2},
    proposed_input::AbstractVector{T},
    options::Dict,
) where {T<:Number}
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]
    replace_invalids = get(options, :replace_invalids, :NoAction)
    dampening = get(options, :dampening, 1.0)
    dampening_with_input = get(options, :dampening_with_input, false)

    # Handle invalid entries (NaN, Inf, missing)
    if eltype(proposed_input) <: Complex
        dodgy_entries =
            (isnan.(real.(proposed_input)) .| isnan.(imag.(proposed_input))) .|
            ismissing.(proposed_input) .|
            (isinf.(real.(proposed_input)) .| isinf.(imag.(proposed_input)))
    else
        dodgy_entries =
            isnan.(proposed_input) .| ismissing.(proposed_input) .| isinf.(proposed_input)
    end

    if sum(dodgy_entries) != 0
        if replace_invalids == :ReplaceElements
            proposed_input[dodgy_entries] = simple_iterate[dodgy_entries]
        elseif replace_invalids == :ReplaceVector
            proposed_input = simple_iterate
        end
    end

    # Apply dampening
    reference_point = dampening_with_input ? inputs[:, end] : simple_iterate
    return @. (dampening * proposed_input) + ((1 - dampening) * reference_point)
end
