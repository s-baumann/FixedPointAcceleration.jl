

# Unified accelerate interface - workspace variants handled by dispatch
accelerate(method::Simple, st::IterationState, cfg::FixedPointConfig, ::Nothing) = st.fx

function accelerate(method::Aitken, st::IterationState, cfg::FixedPointConfig, ::Nothing)
    # Need at least 2 previous (x, fx) pairs => 3 total evaluations
    if length(st.history_x) < 3
        return st.fx
    end
    f1 = st.history_fx[end - 2]
    f2 = st.history_fx[end - 1]
    f3 = st.history_fx[end]
    # Apply elementwise Aitken Δ²: x - (Δx)^2 / (Δ² x)
    Δ1 = f2 .- f1
    Δ2 = f3 .- f2
    denom = Δ2 .- Δ1
    # Avoid division by zero / near-zero -> fallback to latest f(x). For complex, use eps of real type.
    T = eltype(denom)
    ϵ = T <: Real ? eps(T) : eps(float(real(T)))
    mask = abs.(denom) .> 10 * ϵ
    accelerated = similar(f3)
    accelerated .= f3
    accelerated[mask] = f1[mask] .- (Δ1[mask] .^ 2) ./ denom[mask]
    return accelerated
end

# function accelerate(method::Anderson, st::IterationState, cfg::FixedPointConfig, ws=nothing)
#     return ws === nothing ? _anderson_fallback(method, st, cfg) : _anderson_workspace(method, st, cfg, ws)
# end

function accelerate(method::Anderson, st::IterationState, cfg::FixedPointConfig, ::Nothing)
    m = method.m
    k = length(st.history_x)
    if k < 2
        return st.fx
    end
    w = min(m, k)
    start_idx = k - w + 1
    rs = [st.history_fx[start_idx + j - 1] .- st.history_x[start_idx + j - 1] for j in 1:w]
    w < 2 && return st.fx
    ΔR = hcat((rs[i + 1] .- rs[i] for i in 1:(w - 1))...)
    F = ΔR
    Q = qr(F)
    γ = Q \ rs[end]
    if any(isnan.(γ)) || any(isinf.(γ))
        return st.fx
    end
    fxs = [st.history_fx[start_idx + j - 1] for j in 1:w]
    ΔF = hcat((fxs[i + 1] .- fxs[i] for i in 1:(w - 1))...)
    proposed = fxs[end] .- (ΔF * γ)
    return proposed
end

function accelerate(
    method::Anderson, st::IterationState, cfg::FixedPointConfig, ws::Workspace{T}
) where {T}
    m = method.m
    k = length(st.history_x)
    if k < 2
        return st.fx
    end
    w = min(m, k)
    start_idx = k - w + 1
    fxs = [st.history_fx[start_idx + j - 1] for j in 1:w]
    n = length(st.x)
    # Fill residuals columns
    @inbounds for j in 1:w
        idx = start_idx + j - 1
        xj = st.history_x[idx]
        fxj = fxs[j]
        for i in 1:n
            ws.residuals[i, j] = fxj[i] - xj[i]
        end
    end
    if w < 2
        return st.fx
    end
    # Build delta matrices
    cols = w - 1
    @inbounds for j in 1:cols
        for i in 1:n
            ws.delta_resids[i, j] = ws.residuals[i, j + 1] - ws.residuals[i, j]
            ws.delta_outputs[i, j] = (fxs[j + 1][i] - fxs[j][i])
        end
    end
    # Least squares solve ΔR * γ = r_k
    r_last = view(ws.residuals, :, w)
    ΔR = @view ws.delta_resids[:, 1:cols]
    begin
        # Use \; fallback to pinv if needed
        try
            ws.coeffs[1:cols] = ΔR \ r_last
        catch
            ws.coeffs[1:cols] = pinv(Matrix(ΔR)) * r_last
        end
    end
    # proposed = f_k - ΔF * γ
    ΔF = @view ws.delta_outputs[:, 1:cols]
    for i in 1:n
        acc = zero(T)
        @inbounds for j in 1:cols
            acc += ΔF[i, j] * ws.coeffs[j]
        end
        ws.proposed[i] = fxs[end][i] - acc
    end
    return ws.proposed
end
