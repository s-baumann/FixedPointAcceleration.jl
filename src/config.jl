Base.@kwdef struct FixedPointConfig{M,RT,TT}
    metric::M = (x, fx) -> maximum(abs.(fx .- x))
    threshold::TT = 1e-10
    max_iters::Int = 1_000
    relaxation::RT = 1.0                     # damping factor ω
    relaxation_reference::Symbol = :Input    # :Input or :Simple (f(x))
    replace_invalids::Symbol = :NoAction     # :NoAction | :ReplaceElements | :ReplaceVector
    report::Bool = false
    sigfigs::Int = 10
    history_window::Int = -1                 # -1 => keep all, 0 => only last, k => last k
    stall_window::Int = 0                    # 0 disables stall detection
    stall_tolerance::Float64 = 1e-12         # relative improvement threshold
    divergence_factor::Float64 = 1e6         # residual factor triggering divergence
    function FixedPointConfig(
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
        max_iters > 0 || throw(ArgumentError("max_iters must be positive"))
        threshold >= 0 || throw(ArgumentError("threshold must be ≥ 0"))
        (relaxation > 0 && relaxation <= 1) ||
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
        stall_tolerance >= 0 || throw(ArgumentError("stall_tolerance must be ≥ 0"))
        divergence_factor > 0 || throw(ArgumentError("divergence_factor must be > 0"))
        new{typeof(metric),typeof(relaxation),typeof(threshold)}(
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
end

default_config() = FixedPointConfig()
