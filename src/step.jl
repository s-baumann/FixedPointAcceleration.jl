
struct ProgressTracker{R<:Number}
    best_rnorm::R
    best_iter::Int
end

"""Internal convergence status evaluation shared by step variants."""
function evaluate_status!(
    st::IterationState,
    cfg::FixedPointConfig,
    method::AbstractAccelerationMethod,
    iter::Int,
    rnorm::Real,
    tracker::ProgressTracker{R},
) where {R<:Number}
    _maybe_report(cfg, method, iter, rnorm)
    if rnorm <= cfg.threshold
        return :Converged, ProgressTracker{R}(convert(R, rnorm), iter)
    end
    if rnorm > cfg.divergence_factor * st.initial_residual_norm
        return :Diverged, tracker
    end
    if cfg.stall_window > 0
        best_rnorm = tracker.best_rnorm
        improved = rnorm < best_rnorm - cfg.stall_tolerance * best_rnorm
        if improved
            tracker = ProgressTracker{R}(convert(R, rnorm), iter)
        elseif iter - tracker.best_iter >= cfg.stall_window
            return :Stalled, tracker
        end
    end
    return :Continue, tracker
end

function compute_simple_next!(st::IterationState{E}) where {E}
    prev_simple = st.history_simple_x[end]
    simple_next = st.callbacks.f_apply_simple!(prev_simple)
    push!(st.history_simple_x, simple_next)
    return nothing
end

function build_proposal(
    method::AbstractAccelerationMethod, st::IterationState, cfg::FixedPointConfig, ws
)
    st.residual = st.fx .- st.x
    proposed = accelerate(method, st, cfg, ws)
    # Check if acceleration was actually applied (proposed != st.fx means acceleration was used)
    do_accel = proposed !== st.fx
    return do_accel, proposed
end

function finalize_iteration!(
    st::IterationState{E},
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    iter::Int,
    x_new::AbstractVector{E},
    do_accel::Bool,
    tracker::ProgressTracker{R},
    rnorm_prev,
) where {E,R<:Number}
    fx_new = st.callbacks.f_apply_current!(x_new)
    if length(fx_new) != length(x_new)
        return :FunctionSizeMismatch, rnorm_prev, tracker
    end
    st.callbacks.finalize_x!(st, x_new, fx_new)
    update_history_window!(method, cfg, st, do_accel)
    rnorm = cfg.metric(st.x, st.fx)
    status, tracker = evaluate_status!(st, cfg, method, iter, rnorm, tracker)
    return status, rnorm, tracker
end

"""Core step implementation parameterized by evaluation callbacks."""
function step!(
    st::IterationState{E},
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    ws,
    iter::Int,
    tracker::ProgressTracker{R},
    rnorm_prev,
) where {E,R<:Number}
    compute_simple_next!(st)
    do_accel, proposed = build_proposal(method, st, cfg, ws)
    x_new = _apply_relaxation(method, proposed, st, cfg, do_accel)
    return finalize_iteration!(st, method, cfg, iter, x_new, do_accel, tracker, rnorm_prev)
end
