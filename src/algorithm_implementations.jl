"""
Internal function to compute the raw proposed input (before post-processing).
This is where the algorithm-specific dispatch happens.
The actual implementations are in the individual algorithm files.
"""
function _compute_proposed_input end

"""
    fixed_point_new_input(inputs, outputs, algorithm, options, simple_start_index)

Compute the next input using the provided algorithm and global `FixedPointOptions`.
Avoids per-iteration `Dict` creation by passing the structured options and the
changing `simple_start_index` explicitly.
"""
function fixed_point_new_input(
    inputs::AbstractArray{T,2},
    outputs::AbstractArray{T,2},
    algorithm::FixedPointAlgorithm,
    options::FixedPointOptions,
    simple_start_index::Int,
) where {T<:Number}
    proposed_input = _compute_proposed_input(
        inputs, outputs, algorithm, options, simple_start_index
    )
    return _apply_post_processing(inputs, outputs, proposed_input, options)
end

"""
Apply post-processing: replacement strategies and dampening.
"""
function _apply_post_processing(
    inputs::AbstractArray{T,2},
    outputs::AbstractArray{T,2},
    proposed_input::AbstractVector{T},
    options::FixedPointOptions,
) where {T<:Number}
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]
    replace_invalids = options.stability.replace_invalids
    dampening = options.stability.dampening
    dampening_with_input = options.stability.dampening_with_input

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
