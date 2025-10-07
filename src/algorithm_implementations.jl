# New dispatch-based interface
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
    options::Dict
) where {T<:Number}
    # Get the raw proposed input from the algorithm-specific implementation
    proposed_input = _compute_proposed_input(inputs, outputs, algorithm, options)

    # Apply dampening and replacement strategies
    return _apply_post_processing(inputs, outputs, proposed_input, options)
end

"""
Internal function to compute the raw proposed input (before post-processing).
This is where the algorithm-specific dispatch happens.
"""
_compute_proposed_input(inputs, outputs, ::Simple, options) = outputs[:, end]

_compute_proposed_input(inputs, outputs, alg::Anderson, options) =
    _anderson_acceleration(inputs, outputs, alg, options)

_compute_proposed_input(inputs, outputs, ::Aitken, options) =
    _aitken_acceleration(inputs, outputs, options)

_compute_proposed_input(inputs, outputs, ::Newton, options) =
    _newton_acceleration(inputs, outputs, options)

_compute_proposed_input(inputs, outputs, alg::Union{MPE,RRE}, options) =
    _polynomial_extrapolation(inputs, outputs, alg, options)

_compute_proposed_input(inputs, outputs, alg::Union{VEA,SEA}, options) =
    _epsilon_extrapolation(inputs, outputs, alg, options)

# Algorithm-specific implementations
function _anderson_acceleration(
    inputs::AbstractArray{T,2},
    outputs::AbstractArray{T,2},
    alg::Anderson,
    options
) where {T<:Number}
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    if completed_iters < 2
        if get(options, :print_reports, false)
            print("                           Used:", lpad(0, 3), " lags. ")
        end
        return simple_iterate
    end

    vector_length = size(outputs)[1]
    M = min(alg.maxM - 1, completed_iters - 1, vector_length)

    recent_outputs = outputs[:, (completed_iters - M):completed_iters]
    recent_inputs = inputs[:, (completed_iters - M):completed_iters]
    Resid = recent_outputs .- recent_inputs
    DeltaOutputs = recent_outputs[:, 2:(M + 1)] .- recent_outputs[:, 1:M]
    DeltaResids = Resid[:, 2:(M + 1)] .- Resid[:, 1:M]
    LastResid = Resid[:, M + 1]
    LastOutput = recent_outputs[:, M + 1]
    Coeffs = repeat([NaN], size(DeltaOutputs)[2])
    ConditionNumber = NaN

    while any(isnan.(Coeffs))
        if isempty(DeltaResids)
            break
        end
        ConditionNumber = cond(DeltaResids)
        if ConditionNumber > alg.condition_threshold
            M = M - 1
            DeltaOutputs = DeltaOutputs[:, 2:(M + 1)]
            DeltaResids = DeltaResids[:, 2:(M + 1)]
            Coeffs = repeat([NaN], size(DeltaOutputs)[2])
            continue
        end
        # Handle complex numbers by using pinv instead of GLM.fit
        if eltype(DeltaResids) <: Complex
            Coeffs = pinv(DeltaResids) * LastResid
        else
            Fit = fit(LinearModel, hcat(DeltaResids), LastResid)
            Coeffs = Fit.pp.beta0
        end
        if any(isnan.(Coeffs))
            M = M - 1
            if (M < 1.5)
                if get(options, :print_reports, false)
                    print("                          Used:", lpad(0, 3), " lags. ")
                end
                break
            end
            DeltaOutputs = DeltaOutputs[:, 2:(M + 1)]
            DeltaResids = DeltaResids[:, 2:(M + 1)]
        end
    end

    if isempty(Coeffs)
        if get(options, :print_reports, false)
            print(
                "Condition number is ",
                lpad("NaN", 5),
                ". Used:",
                lpad(0, 3),
                " lags. ",
            )
        end
        return repeat([NaN], vector_length)
    else
        if get(options, :print_reports, false)
            print(
                "Condition number is ",
                lpad(round(ConditionNumber; sigdigits=2), 5),
                ". Used:",
                lpad(M + 1, 3),
                " lags. ",
            )
        end
        return LastOutput .- vec(DeltaOutputs * Coeffs)
    end
end

function _aitken_acceleration(inputs, outputs, options)
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

function _newton_acceleration(inputs, outputs, options)
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

function _polynomial_extrapolation(inputs, outputs, alg, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    simple_iterates_matrix = put_together_without_jumps(inputs, outputs)
    if (size(simple_iterates_matrix)[2] % alg.extrapolation_period == 0)
        return PolynomialExtrapolation(simple_iterates_matrix, algorithm_to_symbol(alg))
    else
        return simple_iterate
    end
end

function _epsilon_extrapolation(inputs, outputs, alg, options)
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    simple_iterates_matrix = put_together_without_jumps(inputs, outputs)
    if (size(simple_iterates_matrix)[2] % alg.extrapolation_period) == 0
        return EpsilonExtrapolation(simple_iterates_matrix, algorithm_to_symbol(alg))
    else
        return simple_iterate
    end
end

"""
Apply post-processing: replacement strategies and dampening.
"""
function _apply_post_processing(
    inputs::AbstractArray{T,2},
    outputs::AbstractArray{T,2},
    proposed_input::AbstractVector{T},
    options::Dict
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

# End of algorithm implementations
