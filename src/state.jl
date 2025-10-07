mutable struct IterationState{T}
    x::Vector{T}
    fx::Vector{T}
    residual::Vector{T}
    history_x::Vector{Vector{T}}
    history_fx::Vector{Vector{T}}
    history_simple_x::Vector{Vector{T}}  # pure Picard iterates (f applied without acceleration)
    iter::Int
    initial_residual_norm::Float64
end
function IterationState{E}(x::AbstractVector{E}, fx::AbstractVector{E}) where {E}
    return IterationState{E}(
        collect(x),
        collect(fx),
        fx .- x,
        [collect(x)],
        [collect(fx)],
        [collect(fx)],
        0,
        maximum(abs.(fx .- x)),
    )
end
function init_state(x0::AbstractVector{T}, f::Function) where {T}
    fx0 = f(x0)
    length(fx0) == length(x0) || throw(ArgumentError("Function output length mismatch"))
    residual0 = fx0 .- x0
    r0 = maximum(abs.(residual0))
    # Simple history starts with initial iterate x0; future pushes add f(x_k)
    return IterationState{T}(
        collect(x0),
        collect(fx0),
        residual0,
        [collect(x0)],
        [collect(fx0)],
        [collect(x0)],
        0,
        r0,
    )
end
