function _iterates_from_history(st::IterationState)
    # Use pure Picard iterates (f(x_{k-1})) history for polynomial methods.
    # These are stored in history_simple_x and exclude accelerated outputs,
    # preserving the theoretical polynomial structure required by MPE/RRE.
    return st.history_simple_x
end

# --- Allocation-heavy fallback versions (used if no workspace provided) ---
function accelerate_poly(method::MPE, st::IterationState, cfg::FixedPointConfig)
    seq = _iterates_from_history(st)
    k = length(seq)
    k < 3 && return st.history_simple_x[end]
    (k % method.period) != 0 && return st.history_simple_x[end]
    iters = seq
    n = length(iters[1])
    mcols = k
    M = Matrix{eltype(st.x)}(undef, n, mcols)
    @inbounds for j in 1:mcols
        M[:, j] = iters[j]
    end
    total_columns = mcols
    old_differences = M[:, 2:(total_columns - 1)] .- M[:, 1:(total_columns - 2)]
    last_difference = M[:, total_columns] .- M[:, total_columns - 1]
    cvec = -(pinv(old_differences) * last_difference)
    cvec = vcat(cvec, one(eltype(cvec)))
    s = sum(cvec)
    return (M[:, 2:total_columns] * cvec) ./ s
end

function accelerate_poly(method::RRE, st::IterationState, cfg::FixedPointConfig)
    seq = _iterates_from_history(st)
    k = length(seq)
    k < 4 && return st.history_simple_x[end]
    (k % method.period) != 0 && return st.history_simple_x[end]
    iters = seq
    n = length(iters[1])
    mcols = k
    # Scalar fallback: RRE is poorly conditioned in 1D; defer to simple iterate
    n == 1 && return st.fx
    M = Matrix{eltype(st.x)}(undef, n, mcols)
    @inbounds for j in 1:mcols
        M[:, j] = iters[j]
    end
    first_column = M[:, 1]
    differences = M[:, 2:mcols] - M[:, 1:(mcols - 1)]             # n x (mcols-1)
    # Legacy: second differences use all adjacent pairs -> size n x (mcols-2)
    second_differences = differences[:, 2:end] - differences[:, 1:(end - 1)]
    first_difference = differences[:, 1]                           # n
    base_diffs = differences[:, 1:(end - 1)]                       # n x (mcols-2)
    # Conditioning guard via SVD
    T = eltype(M)
    begin
        sv = try
            svd(second_differences; full=false)
        catch
            if get(ENV, "FPA_RRE_DEBUG", "") == "1"
                println("[RRE2 alloc] iteration=$(mcols) svd failure -> fallback")
            end
            return st.fx
        end
        smax = maximum(sv.S)
        smin = minimum(sv.S)
        if smin == zero(smin) || smax == zero(smax)
            if get(ENV, "FPA_RRE_DEBUG", "") == "1"
                println("[RRE2 alloc] iteration=$(mcols) zero singular value -> fallback")
            end
            return st.fx
        end
        condnum = smax / smin
        if condnum > 1e12
            if get(ENV, "FPA_RRE_DEBUG", "") == "1"
                println(
                    "[RRE2 alloc] iteration=$(mcols) cond=$(condnum) > 1e12 -> fallback"
                )
            end
            return st.fx
        end
        inv_second = sv.V * Diagonal(one(T) ./ sv.S) * sv.U'
        proposed = first_column - (base_diffs * (inv_second * first_difference))
        if get(ENV, "FPA_RRE_DEBUG", "") == "1"
            println(
                "[RRE2 alloc] iteration=$(mcols) cond=$(condnum) ||Δ||=$(maximum(abs.(first_difference))) prop_norm=$(maximum(abs.(proposed)))",
            )
        end
    end
    if any(!isfinite, proposed)
        return st.fx
    end
    return proposed
end

# --- Workspace-based versions (reuse matrices, reduce allocations) ---

function _ensure_capacity!(ws::Workspace{T}, n::Int, k::Int) where {T}
    if size(ws.residuals, 1) != n || size(ws.residuals, 2) < k
        newcols = max(k, size(ws.residuals, 2))
        ws.residuals = zeros(T, n, newcols)
    end
    if size(ws.delta_resids, 1) != n || size(ws.delta_resids, 2) < max(k - 2, 1)
        ws.delta_resids = zeros(T, n, max(k - 2, 1))
    end
    if size(ws.delta_outputs, 1) != n || size(ws.delta_outputs, 2) < max(k - 3, 1)
        ws.delta_outputs = zeros(T, n, max(k - 3, 1))
    end
    if length(ws.coeffs) < max(k - 1, 1)
        ws.coeffs = zeros(T, max(k - 1, 1))
    end
    if length(ws.proposed) != n
        ws.proposed = zeros(T, n)
    end
end

function accelerate_poly(
    method::MPE, st::IterationState, cfg::FixedPointConfig, ws::Workspace{T}
) where {T}
    seq = _iterates_from_history(st)
    k = length(seq)
    k < 3 && return st.history_simple_x[end]
    (k % method.period) != 0 && return st.history_simple_x[end]
    n = length(seq[1])
    _ensure_capacity!(ws, n, k)
    @inbounds for j in 1:k
        col = seq[j]
        for i in 1:n
            ws.residuals[i, j] = col[i]
        end
    end
    # old differences into delta_resids (n x (k-2))
    cols_old = k - 2
    @inbounds for j in 1:cols_old
        for i in 1:n
            ws.delta_resids[i, j] = ws.residuals[i, j + 1] - ws.residuals[i, j]
        end
    end
    # last difference into proposed (temp)
    @inbounds for i in 1:n
        ws.proposed[i] = ws.residuals[i, k] - ws.residuals[i, k - 1]
    end
    old_differences = @view ws.delta_resids[:, 1:cols_old]
    last_difference = ws.proposed
    cvec = try
        -(pinv(old_differences) * last_difference)
    catch
        # Fallback if pinv fails unexpectedly
        return st.fx
    end
    # reuse coeffs vector to append 1 (allocate once if length mismatch)
    need = length(cvec) + 1
    if length(ws.coeffs) < need
        ws.coeffs = zeros(T, need)
    end
    @inbounds for i in 1:length(cvec)
        ws.coeffs[i] = cvec[i]
    end
    ws.coeffs[need] = one(T)
    s = zero(T)
    @inbounds for i in 1:need
        s += ws.coeffs[i]
    end
    # proposed final result: (M[:,2:k] * coeffs)/s
    @inbounds for i in 1:n
        acc = zero(T)
        for j in 2:k
            acc += ws.residuals[i, j] * ws.coeffs[j - 1]
        end
        ws.proposed[i] = acc / s
    end
    return ws.proposed
end

function accelerate_poly(
    method::RRE, st::IterationState, cfg::FixedPointConfig, ws::Workspace{T}
) where {T}
    seq = _iterates_from_history(st)
    k = length(seq)
    k < 4 && return st.history_simple_x[end]
    (k % method.period) != 0 && return st.history_simple_x[end]
    n = length(seq[1])
    n == 1 && return st.fx
    _ensure_capacity!(ws, n, k)
    @inbounds for j in 1:k
        col = seq[j]
        for i in 1:n
            ws.residuals[i, j] = col[i]
        end
    end
    # differences into delta_resids (n x (k-1)) -- grow if needed
    if size(ws.delta_resids, 2) < k - 1
        ws.delta_resids = zeros(T, n, k - 1)
    end
    @inbounds for j in 1:(k - 1)
        for i in 1:n
            ws.delta_resids[i, j] = ws.residuals[i, j + 1] - ws.residuals[i, j]
        end
    end
    # second differences into delta_outputs (n x (k-2)) adjacent pairs
    if size(ws.delta_outputs, 2) < k - 2
        ws.delta_outputs = zeros(T, n, k - 2)
    end
    @inbounds for j in 1:(k - 2)
        for i in 1:n
            ws.delta_outputs[i, j] = ws.delta_resids[i, j + 1] - ws.delta_resids[i, j]
        end
    end
    first_column = view(ws.residuals, :, 1)
    differences = @view ws.delta_resids[:, 1:(k - 1)]          # n x (k-1)
    second_differences = @view ws.delta_outputs[:, 1:(k - 2)]  # n x (k-2)
    first_difference = view(ws.delta_resids, :, 1)             # n
    base_diffs = @view ws.delta_resids[:, 1:(k - 2)]           # n x (k-2)
    # Conditioning guard via SVD
    begin
        sv = try
            svd(second_differences; full=false)
        catch
            if get(ENV, "FPA_RRE_DEBUG", "") == "1"
                println("[RRE2 ws] k=$(k) svd failure -> fallback")
            end
            return st.fx
        end
        smax = maximum(sv.S)
        smin = minimum(sv.S)
        if smin == zero(smin) || smax == zero(smax)
            if get(ENV, "FPA_RRE_DEBUG", "") == "1"
                println("[RRE2 ws] k=$(k) zero singular value -> fallback")
            end
            return st.fx
        end
        condnum = smax / smin
        if condnum > 1e12
            if get(ENV, "FPA_RRE_DEBUG", "") == "1"
                println("[RRE2 ws] k=$(k) cond=$(condnum) > 1e12 -> fallback")
            end
            return st.fx
        end
        inv_second = sv.V * Diagonal(one(T) ./ sv.S) * sv.U'
        tmp = inv_second * first_difference
        if get(ENV, "FPA_RRE_DEBUG", "") == "1"
            println(
                "[RRE2 ws] k=$(k) cond=$(condnum) ||Δ||=$(maximum(abs.(first_difference))) tmp_norm=$(maximum(abs.(tmp)))",
            )
        end
    end
    # proposed = first_column - base_diffs * tmp
    @inbounds for i in 1:n
        acc = zero(T)
        for j in 1:(k - 2)
            acc += base_diffs[i, j] * tmp[j]
        end
        ws.proposed[i] = first_column[i] - acc
    end
    if any(!isfinite, ws.proposed)
        return st.fx
    end
    return ws.proposed
end

# Unified accelerate interface for polynomial methods with built-in acceleration decision
function accelerate(method::Union{RRE,MPE}, st::IterationState, cfg::FixedPointConfig, ws=nothing)
    # Check if we should accelerate (embedded logic from _should_accelerate)
    ksimple = length(st.history_simple_x)
    should_accel = (ksimple % method.period) == 0 && ksimple >= _min_history(method)

    if !should_accel
        return st.fx  # Return st.fx to indicate no acceleration
    end

    # Proceed with acceleration
    return ws === nothing ? accelerate_poly(method, st, cfg) : accelerate_poly(method, st, cfg, ws)
end

# Trait system for minimum history requirements
_min_history(::MPE) = 3
_min_history(::RRE) = 4

# Polynomial methods don't use relaxation - they return the proposed value directly
