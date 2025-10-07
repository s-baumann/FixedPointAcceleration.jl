mutable struct Workspace{T}
    residuals::Matrix{T}
    delta_resids::Matrix{T}
    delta_outputs::Matrix{T}
    coeffs::Vector{T}
    proposed::Vector{T}
end

function init_workspace(::Type{T}, n::Int, m::Int) where {T}
    m < 1 && (m = 1)
    residuals = zeros(T, n, m)
    k = max(m - 1, 1)
    delta_resids = zeros(T, n, k)
    delta_outputs = zeros(T, n, k)
    coeffs = zeros(T, k)
    proposed = similar(residuals[:, 1])
    return Workspace(residuals, delta_resids, delta_outputs, coeffs, proposed)
end

function _maybe_workspace(
    ::Type{T}, st::IterationState{T}, method::AbstractAccelerationMethod
) where {T}
    if method isa Anderson
        return init_workspace(T, length(st.x), method.m)
    elseif method isa MPE || method isa RRE || method isa VEA || method isa SEA
        # polynomial methods need dynamic width; start with small workspace sized to current history
        return init_workspace(T, length(st.x), max(length(st.history_x), 3))
    else
        return nothing
    end
end

function _enforce_history_window!(st::IterationState, cfg::FixedPointConfig)
    hw = cfg.history_window
    hw == -1 && return nothing
    if hw == 0
        st.history_x = [st.history_x[end]]
        st.history_fx = [st.history_fx[end]]
        st.history_simple_x = [st.history_simple_x[end]]  # Ensure simple history is maintained
    else
        lenx = length(st.history_x)
        if lenx > hw
            start = lenx - hw + 1
            start = start < 1 ? 1 : start
            st.history_x = st.history_x[start:lenx]
            st.history_fx = st.history_fx[start:lenx]
        end
        lens = length(st.history_simple_x)
        if lens > hw
            sstart = lens - hw + 1
            sstart = sstart < 1 ? 1 : sstart
            st.history_simple_x = st.history_simple_x[sstart:lens]
        end
    end
end

# Synchronize simple history after a successful polynomial (MPE/RRE) acceleration
_sync_poly_history!(::AbstractAccelerationMethod, st::IterationState, do_accel::Bool) = nothing
function _sync_poly_history!(::Union{MPE,RRE}, st::IterationState, do_accel::Bool)
    do_accel || return nothing
    st.history_simple_x[end] = copy(st.x)
    st.history_simple_x = [st.history_simple_x[end]]
end
