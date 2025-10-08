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

# Synchronize simple history after a successful polynomial (MPE/RRE) acceleration
function sync_poly_history!(
    ::AbstractAccelerationMethod, st::IterationState, cfg::FixedPointConfig, do_accel::Bool
)
    return nothing
end

function sync_poly_history!(
    ::Union{MPE,RRE}, st::IterationState, cfg::FixedPointConfig, do_accel::Bool
)
    do_accel || return nothing
    empty!(st.history_simple_x)
    push!(st.history_simple_x, st.x)
    return nothing
end

function push_history!(st::IterationState, cfg::FixedPointConfig)
    push!(st.history_x, st.x)
    push!(st.history_fx, st.fx)
    return nothing
end

function update_history_window!(
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    st::IterationState,
    do_accel::Bool,
)
    sync_poly_history!(method, st, cfg, do_accel)
    push_history!(st, cfg)
    return nothing
end
