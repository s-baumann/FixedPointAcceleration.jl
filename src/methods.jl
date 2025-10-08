abstract type AbstractAccelerationMethod end

struct Simple <: AbstractAccelerationMethod end

struct Anderson{T<:Real} <: AbstractAccelerationMethod
    m::Int
    λ::T
    function Anderson(; m::Int=5, λ::T=1e-10) where {T<:Real}
        m >= 1 || throw(ArgumentError("m must be ≥ 1"))
        λ >= 0 || throw(ArgumentError("λ must be ≥ 0"))
        new{T}(m, λ)
    end
end

"""Aitken Δ² method (element-wise)."""
struct Aitken <: AbstractAccelerationMethod end

"""Minimal Polynomial Extrapolation (core2)."""
struct MPE <: AbstractAccelerationMethod
    period::Int
    function MPE(; period::Int=7)
        period >= 1 || throw(ArgumentError("period must be ≥ 1"))
        new(period)
    end
end

"""Reduced Rank Extrapolation (core2)."""
struct RRE <: AbstractAccelerationMethod
    period::Int
    function RRE(; period::Int=7)
        period >= 1 || throw(ArgumentError("period must be ≥ 1"))
        new(period)
    end
end

"""Vector Epsilon Algorithm (core2)."""
struct VEA <: AbstractAccelerationMethod
    period::Int
    function VEA(; period::Int=7)
        period >= 1 || throw(ArgumentError("period must be ≥ 1"))
        new(period)
    end
end

"""Scalar Epsilon Algorithm (core2)."""
struct SEA <: AbstractAccelerationMethod
    period::Int
    function SEA(; period::Int=7)
        period >= 1 || throw(ArgumentError("period must be ≥ 1"))
        new(period)
    end
end

# Trait system for relaxation behavior
abstract type RelaxationTrait end
struct UsesRelaxation <: RelaxationTrait end
struct NoRelaxation <: RelaxationTrait end

# Default: most methods use relaxation
_relaxation_trait(::AbstractAccelerationMethod) = UsesRelaxation()
# Polynomial methods don't use relaxation
_relaxation_trait(::Union{MPE,RRE}) = NoRelaxation()

# Trait for methods that use periods
abstract type PeriodTrait end
struct HasPeriod <: PeriodTrait end
struct NoPeriod <: PeriodTrait end

# Default: no period
_period_trait(::AbstractAccelerationMethod) = NoPeriod()
# Methods with periods
_period_trait(::Union{MPE,RRE,VEA,SEA}) = HasPeriod()
