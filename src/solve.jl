struct FixedPointSolution{E,R,M}
    fixed_point::Vector{E}
    residual_norm::R
    iterations::Int
    status::Symbol
    method::M
    history::Vector{Vector{E}}
end

# Apply relaxation / blending rules to proposed iterate.

function solve(
    f::Function,
    x0::AbstractVector{E};
    method::AbstractAccelerationMethod=Simple(),
    cfg::FixedPointConfig=FixedPointConfig(),
) where {E}
    st = init_state(x0, f)
    ws = init_workspace(E, st, method)

    rnorm = cfg.metric(st.x, st.fx)
    _maybe_report(cfg, method, 0, rnorm)
    if rnorm <= cfg.threshold
        return FixedPointSolution(st.x, rnorm, 0, :Converged, method, st.history_x)
    end

    return _run_iterations!(st, method, cfg, ws, rnorm)
end


"""In-place variant expecting `f!(out, x)` writing f(x) to out."""
function solve!(
    f!::Function,
    x::AbstractVector{E};
    method::AbstractAccelerationMethod=Simple(),
    cfg::FixedPointConfig=FixedPointConfig(),
) where {E}
    fx = similar(x)
    f!(fx, x)
    st = init_state!(x, fx, f!)

    ws = init_workspace(E, st, method)
    rnorm = cfg.metric(st.x, st.fx)

    _maybe_report(cfg, method, 0, rnorm)

    if rnorm <= cfg.threshold
        return FixedPointSolution(st.x, rnorm, 0, :Converged, method, st.history_x)
    end

    return _run_iterations!(st, method, cfg, ws, rnorm)
end


function _run_iterations!(
    st::IterationState,
    method::AbstractAccelerationMethod,
    cfg::FixedPointConfig,
    ws,
    rnorm::Real,
)
    tracker = ProgressTracker{typeof(rnorm)}(rnorm, 0)
    for iter in 1:cfg.max_iters
        status, new_rnorm, tracker = step!(st, method, cfg, ws, iter, tracker, rnorm)

        if status === :FunctionSizeMismatch
            return FixedPointSolution(st.x, rnorm, iter - 1, status, method, st.history_x)
        elseif status === :Continue
            rnorm = new_rnorm
        else
            return FixedPointSolution(st.x, new_rnorm, iter, status, method, st.history_x)
        end
    end
    return FixedPointSolution(st.x, rnorm, cfg.max_iters, :MaxIters, method, st.history_x)
end


function _maybe_report(cfg::FixedPointConfig, method, iter::Int, rnorm)
    if cfg.report
        println(
            "Algorithm: ",
            typeof(method).name.name,
            " Iter: ",
            lpad(iter, 5),
            " Residual: ",
            rpad(round(rnorm; sigdigits=cfg.sigfigs), cfg.sigfigs + 6),
            " Time: ",
            now(),
        )
    end
end


function _apply_relaxation(
    method::AbstractAccelerationMethod,
    proposed,
    st::IterationState,
    cfg::FixedPointConfig,
    do_accel::Bool,
)
    if method isa Simple || !do_accel
        return proposed
    end
    ω = cfg.relaxation
    ref = cfg.relaxation_reference === :Input ? st.x : st.fx
    return @. ref + ω * (proposed - ref)
end
