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

"""Return nothing by default for methods that do not need a workspace."""
function init_workspace(
    ::Type{T}, st::IterationState{T}, ::AbstractAccelerationMethod
) where {T}
    return nothing
end

function init_workspace(::Type{T}, st::IterationState{T}, method::Anderson) where {T}
    return init_workspace(T, length(st.x), method.m)
end

function init_workspace(
    ::Type{T}, st::IterationState{T}, ::Union{MPE,RRE,VEA,SEA}
) where {T}
    # polynomial methods need dynamic width; start with small workspace sized to current history
    return init_workspace(T, length(st.x), max(length(st.history_x), 3))
end

function enforce_history_window!(st::IterationState, cfg::FixedPointConfig)
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
function sync_poly_history!(
    ::AbstractAccelerationMethod, st::IterationState, do_accel::Bool
)
    nothing
end
function sync_poly_history!(::Union{MPE,RRE}, st::IterationState, do_accel::Bool)
    do_accel || return nothing
    st.history_simple_x[end] = copy(st.x)
    st.history_simple_x = [st.history_simple_x[end]]
end

function update_history_window!(
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    st::IterationState,
    do_accel::Bool,
)
    sync_poly_history!(method, st, do_accel)
    push_history!(st)
    enforce_history_window!(st, cfg)
    return nothing
end

push_history!(st_local::IterationState) = begin
    push!(st_local.history_x, st_local.x)
    push!(st_local.history_fx, st_local.fx)
end
