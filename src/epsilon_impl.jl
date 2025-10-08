
"""Compute epsilon accelerated proposal (fallback allocation version)."""
function accelerate(method::VEA, st::IterationState, cfg::FixedPointConfig, ws::Nothing)
    return _epsilon_accelerate(method.period, st; vector=true)
end
function accelerate(method::SEA, st::IterationState, cfg::FixedPointConfig, ws::Nothing)
    return _epsilon_accelerate(method.period, st; vector=false)
end

"""Workspace variants: reuse workspace matrix to avoid rebuilding iterate matrix allocations."""
function accelerate(
    method::VEA, st::IterationState, cfg::FixedPointConfig, ws::Workspace
)
    return _epsilon_accelerate_ws(method.period, st, ws; vector=true)
end
function accelerate(
    method::SEA, st::IterationState, cfg::FixedPointConfig, ws::Workspace
)
    return _epsilon_accelerate_ws(method.period, st, ws; vector=false)
end

function _epsilon_accelerate(period::Int, st::IterationState; vector::Bool)
    # Need at least one completed iteration beyond initial to consider.
    length(st.history_x) < 2 && return st.fx
    # Build simple iterate matrix like legacy: [x history..., last fx].
    ncols_inputs = length(st.history_x)
    n = length(st.x)
    iter_cols = ncols_inputs + 1
    if (iter_cols % period) != 0
        return st.fx
    end
    M = Matrix{eltype(st.x)}(undef, n, iter_cols)
    @inbounds for j in 1:ncols_inputs
        M[:, j] = st.history_x[j]
    end
    M[:, iter_cols] = st.fx
    return _epsilon_transform(M; vector=vector)
end

function _epsilon_transform(iterates::AbstractMatrix{T}; vector::Bool) where {T}
    total_cols = size(iterates, 2)
    total_cols == 1 && return iterates[:, 1]
    # Ensure odd number of columns as in legacy implementation.
    if iseven(total_cols)
        iterates = iterates[:, 2:end]
        total_cols = size(iterates, 2)
        total_cols == 1 && return iterates[:, 1]
    end
    mat = iterates
    rows = size(mat, 1)
    prev = zeros(T, rows, total_cols - 1)
    for mc in total_cols:-1:2
        diff = mat[:, 2:mc] .- mat[:, 1:(mc - 1)]
        new = prev + (vector ? _vea_vector_of_inverses(diff) : (1 ./ diff))
        prev = mat[:, 2:(mc - 1)]
        mat = new
    end
    # Retry with fewer columns if NaNs produced (recursive degradation like legacy).
    if any(isnan.(mat))
        total_cols > 3 && return _epsilon_transform(iterates[:, 3:end]; vector=vector)
    end
    return mat[:, 1]
end

function _vea_vector_of_inverses(diff::AbstractMatrix{T}) where {T}
    # Vector epsilon: Moore-Penrose pseudoinverse of each column (treated as vector), transposed.
    if size(diff, 1) == 1
        return 1 ./ diff
    end
    cols = size(diff, 2)
    out = Matrix{T}(undef, size(diff, 1), cols)
    @inbounds for j in 1:cols
        # pinv on a column vector returns a row vector; transpose to column orientation needed here
        out[:, j] = (pinv(diff[:, j]))'  # (1 x n) -> assign as column via ' producing (n x 1)
    end
    return out
end

function _epsilon_accelerate_ws(
    period::Int, st::IterationState, ws::Workspace; vector::Bool
)
    length(st.history_x) < 2 && return st.fx
    ncols_inputs = length(st.history_x)
    n = length(st.x)
    total_cols = ncols_inputs + 1
    (total_cols % period) != 0 && return st.fx
    # Ensure workspace matrix capacity
    if size(ws.residuals, 1) != n || size(ws.residuals, 2) < total_cols
        ws.residuals = zeros(eltype(st.x), n, total_cols)
    end
    @inbounds for j in 1:ncols_inputs
        col = st.history_x[j]
        for i in 1:n
            ws.residuals[i, j] = col[i]
        end
    end
    @inbounds for i in 1:n
        ws.residuals[i, total_cols] = st.fx[i]
    end
    # Use a view; transformation allocates internally but avoids rebuilding M
    Mview = view(ws.residuals, :, 1:total_cols)
    return _epsilon_transform(Mview; vector=vector)
end
