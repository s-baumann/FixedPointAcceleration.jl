struct IterationCallbacks{F1,F2,F3}
    f_apply_current!::F1
    f_apply_simple!::F2
    finalize_x!::F3
end

function history_capacity(cfg::FixedPointConfig)
    cfg.history_window == -1 ? cfg.max_iters + 1 : max(cfg.history_window, 1)
end

function _init_history_buffer(cap::Int, v::Vector{T}) where {T}
    buf = CircularBuffer{Vector{T}}(cap)
    push!(buf, copy(v))
    return buf
end

mutable struct IterationState{T,R<:Real}
    x::Vector{T}
    fx::Vector{T}
    residual::Vector{T}
    history_x::CircularBuffer{Vector{T}}
    history_fx::CircularBuffer{Vector{T}}
    history_simple_x::CircularBuffer{Vector{T}} # pure Picard iterates (f applied without acceleration)
    iter::Int
    initial_residual_norm::R
    callbacks::IterationCallbacks
end

function init_state(x0::AbstractVector{T}, f::Function, history_cap::Int) where {T}
    fx0 = f(x0)
    length(fx0) == length(x0) || throw(ArgumentError("Function output length mismatch"))
    x_copy = collect(x0)
    fx_copy = collect(fx0)
    residual = fx_copy .- x_copy
    r0 = maximum(abs.(residual))
    f_apply_current!(x::AbstractVector{T}) = f(x)
    f_apply_simple!(prev_simple::AbstractVector{T}) = f(prev_simple)
    finalize_x!(
        st_local::IterationState{T}, x_new::AbstractVector{T}, fx_new::AbstractVector{T}
    ) = begin
        st_local.x = x_new
        st_local.fx = fx_new
        st_local.residual = fx_new .- x_new
        nothing
    end
    cb = IterationCallbacks(f_apply_current!, f_apply_simple!, finalize_x!)

    st = IterationState{T,typeof(r0)}(
        x_copy,
        fx_copy,
        residual,
        _init_history_buffer(history_cap, x_copy),
        _init_history_buffer(history_cap, fx_copy),
        _init_history_buffer(history_cap, x_copy),
        0,
        r0,
        cb,
    )
    return st
end

function init_state!(
    x::AbstractVector{T}, fx::AbstractVector{T}, f!::Function, history_cap::Int
) where {T}
    f_apply_current!(x_candidate::AbstractVector{T}) = begin
        f!(fx, x_candidate)
        fx
    end
    f_apply_simple!(prev_simple::AbstractVector{T}) = begin
        tmp = copy(prev_simple)
        f!(tmp, prev_simple)
        tmp
    end
    finalize_x!(
        st_local::IterationState{T}, x_new::AbstractVector{T}, fx_new::AbstractVector{T}
    ) = begin
        x .= x_new
        st_local.x .= x
        st_local.fx .= fx_new
        st_local.residual .= fx_new .- x
        nothing
    end
    cb = IterationCallbacks(f_apply_current!, f_apply_simple!, finalize_x!)

    x_copy = collect(x)
    fx_copy = collect(fx)
    residual = fx_copy .- x_copy
    r0 = maximum(abs.(residual))
    return IterationState{T,typeof(r0)}(
        x_copy,
        fx_copy,
        residual,
        _init_history_buffer(history_cap, x_copy),
        _init_history_buffer(history_cap, fx_copy),
        _init_history_buffer(history_cap, x_copy),
        0,
        r0,
        cb,
    )
    return st
end
