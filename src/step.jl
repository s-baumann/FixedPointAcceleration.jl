
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

function compute_simple_next!(
    st::IterationState{E}
) where {E}
    prev_simple = st.history_simple_x[end]
    simple_next = st.callbacks.f_apply_simple!(prev_simple)
    push!(st.history_simple_x, simple_next)
    return nothing
end

function build_proposal(
    method::AbstractAccelerationMethod,
    st::IterationState,
    cfg::FixedPointConfig,
    ws,
)
    st.residual = st.fx .- st.x
    do_accel = _should_accelerate(method, st)
    proposed = do_accel ? _accelerated_proposal(method, st, cfg, ws) : st.fx
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
    return finalize_iteration!(
        st,
        method,
        cfg,
        iter,
        x_new,
        do_accel,
        tracker,
        rnorm_prev,
    )
end

function update_history_window!(
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    st::IterationState,
    do_accel::Bool,
)
    _sync_poly_history!(method, st, do_accel)
    push_history!(st)
    _enforce_history_window!(st, cfg)
    return nothing
end

push_history!(st_local::IterationState) = begin
    push!(st_local.history_x, st_local.x)
    push!(st_local.history_fx, st_local.fx)
end
