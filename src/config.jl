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

function FixedPointConfig(;
    metric::M = (x, fx) -> maximum(abs.(fx .- x)),
    threshold::TT  = 1e-10,
    max_iters::Integer = 1_000,
    relaxation::RT  = 1.0,
    relaxation_reference::Symbol = :Input,
    replace_invalids::Symbol= :NoAction,
    report::Bool = false,
    sigfigs::Int = 10,
    history_window::Int = -1,
    stall_window::Int = 0,
    stall_tolerance::Real = 1e-12,
    divergence_factor::Real = 1e6,
) where {M,RT<:Real,TT<:Real}
    Int(max_iters) > 0 || throw(ArgumentError("max_iters must be positive"))
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
    Float64(stall_tolerance) >= 0 || throw(ArgumentError("stall_tolerance must be ≥ 0"))
    Float64(divergence_factor) > 0 || throw(ArgumentError("divergence_factor must be > 0"))

    return FixedPointConfig{M,RT,TT}(
        metric,
        threshold,
        Int(max_iters),
        relaxation,
        relaxation_reference,
        replace_invalids,
        report,
        sigfigs,
        history_window,
        stall_window,
        Float64(stall_tolerance),
        Float64(divergence_factor),
    )
end

default_config() = FixedPointConfig()
