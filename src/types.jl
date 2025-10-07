"""
 This contains the results from a function evaluation in the FixedPointAcceleration framework.
 It containing the following fields:
 * `Input_` - The input
 * `Output_` - The output of the  function. May be `missing` if function could not complete without error.
 * `Error_` - A symbol representing what error occured.
"""
struct FunctionEvaluationResult{T<:Number,R}
    input::Vector{T}
    output::Union{Missing,Vector{Missing},Vector{Union{Missing,R}},Vector{R}}
    other_output::Union{Missing,NamedTuple}
    error::Symbol
    function FunctionEvaluationResult(
        input_vec::Vector{T},
        output_val::Missing,
        error_sym::Symbol,
        other_output_val::Union{Missing,NamedTuple}=missing,
    ) where {T<:Number}
        return new{T,T}(input_vec, output_val, other_output_val, error_sym)
    end
    function FunctionEvaluationResult(
        input_vec::Vector{T},
        output_vec::Vector{R},
        error_sym::Symbol,
        other_output_val::Union{Missing,NamedTuple}=missing,
    ) where {T<:Number} where {R<:Union{Missing,<:Number}}
        if R === Missing
            return new{T,T}(input_vec, output_vec, other_output_val, error_sym)
        elseif R <: Number
            return new{T,R}(
                input_vec,
                convert(Vector{Union{Missing,R}}, output_vec),
                other_output_val,
                error_sym,
            )
        else
            return new{T,nonmissingtype(R)}(
                input_vec, output_vec, other_output_val, error_sym
            )
        end
    end
end

"""
 This contains the results from a the fixed_point function.

 It containing the following fields:
 * `FixedPoint_` - The `Vector` with the fixed point that has been found.
 * `Other_Output_` - The other output of the fixedpoint function.
 * `Convergence_` - A real number showing how close the `FixedPoint_` is to the input (to the function) that created it as output.
  * `TerminationCondition_' - Why did the fixedpoint acceleration stop.
  * `Iterations_' - How many iterations were undertaken
  * `ConvergenceVector_` - What is the convergence value at each iteration.
  * `FailedEvaluation_` - Why did the fixedpoint iteration fail (missing if it did not fail)
  * `Inputs_` - What were all of the inputs tried
  * `Outputs_` - What were all the corresponding outputs.
"""
struct FixedPointResults{R<:Number}
    fixed_point::Union{Missing,Vector{R}}
    other_output::Union{Missing,NamedTuple}
    convergence::Union{Missing,Real}
    termination_condition::Symbol
    iterations::Integer
    convergence_vector::Union{Missing,Vector{<:Real}}
    failed_evaluation::Union{Missing,FunctionEvaluationResult}
    inputs::Matrix{R}
    outputs::Matrix{R}
    function FixedPointResults(
        inputs_mat::Matrix{R},
        outputs_mat::Matrix{R},
        termination_condition::Symbol;
        convergence_vector::Union{Missing,Vector{<:Real}}=missing,
        failed_evaluation::Union{Missing,FunctionEvaluationResult}=missing,
        other_output_val::Union{Missing,NamedTuple}=missing,
    ) where {R<:Number}
        num_iterations = size(outputs_mat, 2)
        fixed_point_val = missing
        convergence_val = missing
        if !ismissing(convergence_vector) && !isempty(convergence_vector)
            convergence_val = convergence_vector[num_iterations]
        end
        if termination_condition == :ReachedConvergenceThreshold
            fixed_point_val = outputs_mat[:, num_iterations]
        end
        return new{R}(
            fixed_point_val,
            other_output_val,
            convergence_val,
            termination_condition,
            num_iterations,
            convergence_vector,
            failed_evaluation,
            inputs_mat,
            outputs_mat,
        )
    end
end

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
