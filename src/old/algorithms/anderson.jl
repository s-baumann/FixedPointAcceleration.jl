"""
Anderson acceleration algorithm implementation.
"""

"""
    Anderson(; maxM=10, condition_threshold=1e3)

Anderson acceleration algorithm.

Anderson acceleration takes a linear combination of previous iterates to accelerate convergence.
Unlike other methods, it can use any previous iterates (not necessarily sequential).

# Parameters
- `maxM::Int`: Maximum number of previous iterates to use (default: 10)
- `condition_threshold::Float64`: Condition number threshold for numerical stability (default: 1e3)

# References
Anderson, D.G. (1965). "Iterative Procedures for Nonlinear Integral Equations."
Journal of the ACM 12(4): 547-60.
"""
struct Anderson <: FixedPointAlgorithm
    maxM::Int
    condition_threshold::Float64

    function Anderson(; maxM::Int=10, condition_threshold::Float64=1e3)
        maxM >= 1 || throw(ArgumentError("maxM must be at least 1"))
        condition_threshold >= 1.0 ||
            throw(ArgumentError("condition_threshold must be at least 1.0"))
        new(maxM, condition_threshold)
    end
end

# Algorithm properties
algorithm_name(::Anderson) = "Anderson"
needs_extrapolation_period(::Anderson) = false
get_extrapolation_period(::Anderson) = 1
is_polynomial_method(::Anderson) = false
is_epsilon_method(::Anderson) = false

# Algorithm implementation
"""
Compute the next input using Anderson acceleration.
"""
function _compute_proposed_input(
    inputs::AbstractArray{T,2},
    outputs::AbstractArray{T,2},
    alg::Anderson,
    options::FixedPointOptions,
    simple_start_index::Int,
) where {T<:Number}
    completed_iters = size(outputs)[2]
    simple_iterate = outputs[:, completed_iters]

    if completed_iters < 2
        if options.print_reports
            print("                           Used:", lpad(0, 3), " lags. ")
        end
        return simple_iterate
    end

    vector_length = size(outputs)[1]
    M = min(alg.maxM - 1, completed_iters - 1, vector_length)

    recent_outputs = outputs[:, (completed_iters - M):completed_iters]
    recent_inputs = inputs[:, (completed_iters - M):completed_iters]
    Resid = recent_outputs .- recent_inputs
    DeltaOutputs = recent_outputs[:, 2:(M + 1)] .- recent_outputs[:, 1:M]
    DeltaResids = Resid[:, 2:(M + 1)] .- Resid[:, 1:M]
    LastResid = Resid[:, M + 1]
    LastOutput = recent_outputs[:, M + 1]
    Coeffs = repeat([NaN], size(DeltaOutputs)[2])
    ConditionNumber = NaN

    while any(isnan.(Coeffs))
        if isempty(DeltaResids)
            break
        end
        ConditionNumber = cond(DeltaResids)
        if ConditionNumber > alg.condition_threshold
            M = M - 1
            DeltaOutputs = DeltaOutputs[:, 2:(M + 1)]
            DeltaResids = DeltaResids[:, 2:(M + 1)]
            Coeffs = repeat([NaN], size(DeltaOutputs)[2])
            continue
        end
        if eltype(DeltaResids) <: Complex
            Coeffs = pinv(DeltaResids) * LastResid
        else
            Fit = fit(LinearModel, hcat(DeltaResids), LastResid)
            Coeffs = Fit.pp.beta0
        end
        if any(isnan.(Coeffs))
            M = M - 1
            if (M < 1.5)
                if options.print_reports
                    print("                          Used:", lpad(0, 3), " lags. ")
                end
                break
            end
            DeltaOutputs = DeltaOutputs[:, 2:(M + 1)]
            DeltaResids = DeltaResids[:, 2:(M + 1)]
        end
    end

    if isempty(Coeffs)
        if options.print_reports
            print("Condition number is ", lpad("NaN", 5), ". Used:", lpad(0, 3), " lags. ")
        end
        return repeat([NaN], vector_length)
    else
        if options.print_reports
            print(
                "Condition number is ",
                lpad(round(ConditionNumber; sigdigits=2), 5),
                ". Used:",
                lpad(M + 1, 3),
                " lags. ",
            )
        end
        return LastOutput .- vec(DeltaOutputs * Coeffs)
    end
end
