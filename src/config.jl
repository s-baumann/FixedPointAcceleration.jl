struct FixedPointConfig{M,RT<:Real,TT<:Real}
    metric::M
    threshold::TT
    max_iters::Int
    relaxation::RT
    relaxation_reference::Symbol
    replace_invalids::Symbol
    report::Bool
    sigfigs::Int
    history_window::Int
    stall_window::Int
    stall_tolerance::Float64
    divergence_factor::Float64
end

function FixedPointConfig(
    metric::M,
    threshold::TT,
    max_iters::Integer,
    relaxation::RT,
    relaxation_reference::Symbol,
    replace_invalids::Symbol,
    report::Bool,
    sigfigs::Int,
    history_window::Int,
    stall_window::Int,
    stall_tolerance::Real,
    divergence_factor::Real,
) where {M,RT<:Real,TT<:Real}
    max_iters_val = Int(max_iters)
    max_iters_val > 0 || throw(ArgumentError("max_iters must be positive"))
    threshold >= zero(TT) || throw(ArgumentError("threshold must be ≥ 0"))
    (relaxation > zero(RT) && relaxation <= one(RT)) ||
        throw(ArgumentError("relaxation ω must satisfy 0 < ω ≤ 1"))
    relaxation_reference in (:Input, :Simple) ||
        throw(ArgumentError("relaxation_reference must be :Input or :Simple"))
    replace_invalids in (:NoAction, :ReplaceElements, :ReplaceVector) || throw(
        ArgumentError(
            "replace_invalids must be :NoAction, :ReplaceElements, or :ReplaceVector",
        ),
    )
    sigfigs > 0 || throw(ArgumentError("sigfigs must be positive"))
    history_window >= -1 || throw(ArgumentError("history_window must be -1 or ≥ 0"))
    stall_window >= 0 || throw(ArgumentError("stall_window must be ≥ 0"))
    stall_tolerance_val = Float64(stall_tolerance)
    stall_tolerance_val >= 0 || throw(ArgumentError("stall_tolerance must be ≥ 0"))
    divergence_factor_val = Float64(divergence_factor)
    divergence_factor_val > 0 || throw(ArgumentError("divergence_factor must be > 0"))
    return FixedPointConfig{M,RT,TT}(
        metric,
        threshold,
        max_iters_val,
        relaxation,
        relaxation_reference,
        replace_invalids,
        report,
        sigfigs,
        history_window,
        stall_window,
        stall_tolerance_val,
        divergence_factor_val,
    )
end

function FixedPointConfig(; kwargs...)
    allowed = (
        :metric,
        :threshold,
        :max_iters,
        :relaxation,
        :relaxation_reference,
        :replace_invalids,
        :report,
        :sigfigs,
        :history_window,
        :stall_window,
        :stall_tolerance,
        :divergence_factor,
    )
    for key in keys(kwargs)
        key in allowed || throw(ArgumentError("unexpected keyword argument: $key"))
    end
    return FixedPointConfig(
        get(kwargs, :metric, (x, fx) -> maximum(abs.(fx .- x))),
        get(kwargs, :threshold, 1e-10),
        get(kwargs, :max_iters, 1_000),
        get(kwargs, :relaxation, 1.0),
        get(kwargs, :relaxation_reference, :Input),
        get(kwargs, :replace_invalids, :NoAction),
        get(kwargs, :report, false),
        get(kwargs, :sigfigs, 10),
        get(kwargs, :history_window, -1),
        get(kwargs, :stall_window, 0),
        get(kwargs, :stall_tolerance, 1e-12),
        get(kwargs, :divergence_factor, 1e6),
    )
end

function FixedPointConfig{M,RT,TT}(
    metric::M,
    threshold::TT,
    max_iters::Integer,
    relaxation::RT,
    relaxation_reference::Symbol,
    replace_invalids::Symbol,
    report::Bool,
    sigfigs::Int,
    history_window::Int,
    stall_window::Int,
    stall_tolerance::Real,
    divergence_factor::Real,
) where {M,RT<:Real,TT<:Real}
    return FixedPointConfig(
        metric,
        threshold,
        max_iters,
        relaxation,
        relaxation_reference,
        replace_invalids,
        report,
        sigfigs,
        history_window,
        stall_window,
        stall_tolerance,
        divergence_factor,
    )
end

function FixedPointConfig{M,RT,TT}(; kwargs...) where {M,RT<:Real,TT<:Real}
    cfg = FixedPointConfig(; kwargs...)
    cfg.metric isa M || throw(ArgumentError("metric must be of type $M"))
    cfg.relaxation isa RT || throw(ArgumentError("relaxation must be of type $RT"))
    cfg.threshold isa TT || throw(ArgumentError("threshold must be of type $TT"))
    return FixedPointConfig(
        cfg.metric,
        cfg.threshold,
        cfg.max_iters,
        cfg.relaxation,
        cfg.relaxation_reference,
        cfg.replace_invalids,
        cfg.report,
        cfg.sigfigs,
        cfg.history_window,
        cfg.stall_window,
        cfg.stall_tolerance,
        cfg.divergence_factor,
    )
end

default_config() = FixedPointConfig()
