
"""Single iteration step for the out-of-place solver.

Returns (status, new_rnorm). A status of :Continue indicates that iteration
should proceed. All other status values correspond to terminal solver states.
"""
function _step!(
    f::Function,
    st::IterationState{E},
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    ws,
    metric::Function,
    iter::Int,
    prev_best::Base.RefValue{R},
    prev_best_iter::Base.RefValue{Int},
    rnorm_prev,
) where {E,R}
    # Extend independent simple (Picard) sequence
    prev_simple = st.history_simple_x[end]
    next_simple = f(prev_simple)
    push!(st.history_simple_x, copy(next_simple))

    # Evaluate current iterate
    fx_pre = f(st.x)
    st.fx = fx_pre
    st.residual = fx_pre .- st.x

    do_accel = _should_accelerate(method, st)
    proposed = if (method isa Simple || !do_accel)
        fx_pre
    else
        _accelerated_proposal(method, st, cfg, ws)
    end
    x_new = _apply_relaxation(method, proposed, st, cfg, do_accel)
    fx_new = f(x_new)
    if length(fx_new) != length(x_new)
        return :FunctionSizeMismatch, rnorm_prev
    end

    st.x = x_new
    st.fx = fx_new
    st.residual = fx_new .- x_new
    _sync_poly_history!(method, st, do_accel)
    push!(st.history_x, copy(x_new))
    push!(st.history_fx, copy(fx_new))
    _enforce_history_window!(st, cfg)

    rnorm = metric(x_new, fx_new)
    _maybe_report(cfg, method, iter, rnorm)
    if rnorm <= cfg.threshold
        return :Converged, rnorm
    end
    if rnorm > cfg.divergence_factor * st.initial_residual_norm
        return :Diverged, rnorm
    end
    if cfg.stall_window > 0
        if rnorm < prev_best[] - cfg.stall_tolerance * prev_best[]
            prev_best[] = rnorm
            prev_best_iter[] = iter
        elseif iter - prev_best_iter[] >= cfg.stall_window
            return :Stalled, rnorm
        end
    end
    return :Continue, rnorm
end

"""Single iteration step for the in-place solver.

Returns (status, new_rnorm). Mirrors `_step!` but operates on user-provided
buffers `xbuf` and `fxbuf`.
"""
function _step_inplace!(
    f!::Function,
    xbuf::AbstractVector{E},
    fxbuf::AbstractVector{E},
    st::IterationState{E},
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    ws,
    iter::Int,
    prev_best::Base.RefValue{R},
    prev_best_iter::Base.RefValue{Int},
    rnorm_prev,
) where {E,R}
    prev_simple = st.history_simple_x[end]
    simple_next = copy(prev_simple)
    f!(simple_next, prev_simple)
    push!(st.history_simple_x, copy(simple_next))

    f!(fxbuf, st.x)
    st.fx .= fxbuf
    st.residual .= fxbuf .- st.x

    do_accel = _should_accelerate(method, st)

    proposed = !do_accel ? st.fx : _accelerated_proposal(method, st, cfg, ws)
    x_new = _apply_relaxation(method, proposed, st, cfg, do_accel)

    @. xbuf = x_new
    f!(fxbuf, xbuf)
    st.x .= xbuf
    st.fx .= fxbuf
    st.residual .= fxbuf .- xbuf

    update_hystory_window!(method, cfg, st, do_accel)

    rnorm = cfg.metric(st.x, st.fx)
    _maybe_report(cfg, method, iter, rnorm)
    if rnorm <= cfg.threshold
        return :Converged, rnorm
    end
    if rnorm > cfg.divergence_factor * st.initial_residual_norm
        return :Diverged, rnorm
    end
    if cfg.stall_window > 0
        if rnorm < prev_best[] - cfg.stall_tolerance * prev_best[]
            prev_best[] = rnorm
            prev_best_iter[] = iter
        elseif iter - prev_best_iter[] >= cfg.stall_window
            return :Stalled, rnorm
        end
    end
    return :Continue, rnorm
end

function update_hystory_window!(
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    st::IterationState,
    do_accel::Bool,
)
    _sync_poly_history!(method, st, do_accel)
    push!(st.history_x, copy(st.x))
    push!(st.history_fx, copy(st.fx))
    return _enforce_history_window!(st, cfg)
end
