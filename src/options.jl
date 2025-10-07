"""
Configuration options for convergence criteria.

# Fields
- `metric::Function`: Function that measures convergence (input, output) -> scalar
- `threshold::Real`: Threshold for convergence (default: 1e-10)
- `max_iterations::Integer`: Maximum number of iterations (default: 1000)
"""
struct ConvergenceOptions
    metric::Function
    threshold::Real
    max_iterations::Integer

    function ConvergenceOptions(;
        metric::Function=(input, output) -> maximum(abs.(output .- input)),
        threshold::Real=1e-10,
        max_iterations::Integer=1000,
    )
        threshold >= 0 || throw(ArgumentError("threshold must be non-negative"))
        max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
        new(metric, threshold, max_iterations)
    end
end

"""
Configuration options for numerical stability and error handling.

# Fields
- `dampening::Real`: Dampening parameter (default: 1.0, no dampening)
- `dampening_with_input::Bool`: Apply dampening to input vs output (default: false)
- `replace_invalids::Symbol`: How to handle NaN/Inf: `:NoAction`, `:ReplaceElements`, `:ReplaceVector`
- `quiet_errors::Bool`: Return partial results on error instead of throwing (default: false)
"""
struct StabilityOptions
    dampening::Real
    dampening_with_input::Bool
    replace_invalids::Symbol
    quiet_errors::Bool

    function StabilityOptions(;
        dampening::Real=1.0,
        dampening_with_input::Bool=false,
        replace_invalids::Symbol=:NoAction,
        quiet_errors::Bool=false,
    )
        0 < dampening <= 1 || throw(ArgumentError("dampening must be in (0, 1]"))
        replace_invalids in [:NoAction, :ReplaceElements, :ReplaceVector] || throw(
            ArgumentError(
                "replace_invalids must be :NoAction, :ReplaceElements, or :ReplaceVector",
            ),
        )
        new(dampening, dampening_with_input, replace_invalids, quiet_errors)
    end
end

"""
Configuration options for progress reporting during iteration.

# Fields
- `print_reports::Bool`: Print iteration progress (default: false)
- `reporting_sig_figs::Integer`: Significant figures for progress reports (default: 10)
"""
struct ReportingOptions
    print_reports::Bool
    reporting_sig_figs::Integer

    function ReportingOptions(; print_reports::Bool=false, reporting_sig_figs::Integer=10)
        reporting_sig_figs > 0 ||
            throw(ArgumentError("reporting_sig_figs must be positive"))
        new(print_reports, reporting_sig_figs)
    end
end

"""
Combined configuration options for fixed point iteration.

# Fields
- `convergence::ConvergenceOptions`: Convergence criteria options
- `stability::StabilityOptions`: Numerical stability options
- `reporting::ReportingOptions`: Progress reporting options

# Constructors
- `FixedPointOptions()`: Default options
- `FixedPointOptions(; kwargs...)`: Specify individual parameters
- `FixedPointOptions(convergence, stability, reporting)`: Specify option groups
"""
struct FixedPointOptions
    convergence::ConvergenceOptions
    stability::StabilityOptions
    reporting::ReportingOptions

    # Constructor with individual parameters (has defaults, can be called with no args)
    function FixedPointOptions(;
        # Convergence options
        metric::Function=(input, output) -> maximum(abs.(output .- input)),
        threshold::Real=1e-10,
        max_iterations::Integer=1000,
        # Stability options
        dampening::Real=1.0,
        dampening_with_input::Bool=false,
        replace_invalids::Symbol=:NoAction,
        quiet_errors::Bool=false,
        # Reporting options
        print_reports::Bool=false,
        reporting_sig_figs::Integer=10,
    )
        convergence = ConvergenceOptions(;
            metric=metric, threshold=threshold, max_iterations=max_iterations
        )
        stability = StabilityOptions(;
            dampening=dampening,
            dampening_with_input=dampening_with_input,
            replace_invalids=replace_invalids,
            quiet_errors=quiet_errors,
        )
        reporting = ReportingOptions(;
            print_reports=print_reports, reporting_sig_figs=reporting_sig_figs
        )
        new(convergence, stability, reporting)
    end

    # Constructor with option groups
    function FixedPointOptions(
        convergence::ConvergenceOptions,
        stability::StabilityOptions,
        reporting::ReportingOptions,
    )
        new(convergence, stability, reporting)
    end
end

# Preset configurations
"""
    default_options()

Default fixed point options suitable for most problems.
"""
default_options() = FixedPointOptions()

"""
    robust_options()

Options with enhanced stability for difficult problems.
- Lower dampening (0.7) for more stable iteration
- Replace invalid elements to avoid NaN propagation
- Quiet errors to return partial results
"""
function robust_options()
    FixedPointOptions(; dampening=0.7, replace_invalids=:ReplaceElements, quiet_errors=true)
end

"""
    fast_options()

Options optimized for fast convergence on well-behaved problems.
- Tighter convergence threshold (1e-12)
- Reduced maximum iterations (500)
- No dampening for fastest convergence
"""
function fast_options()
    FixedPointOptions(; threshold=1e-12, max_iterations=500, dampening=1.0)
end

"""
    debug_options()

Options with detailed reporting for debugging convergence issues.
- Print iteration reports
- Higher precision reporting (12 significant figures)
- Replace invalid elements to avoid crashes
"""
function debug_options()
    FixedPointOptions(;
        print_reports=true, reporting_sig_figs=12, replace_invalids=:ReplaceElements
    )
end
