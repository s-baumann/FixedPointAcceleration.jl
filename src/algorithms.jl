"""
Algorithm types for fixed point acceleration methods.

This module defines the type hierarchy for different fixed point acceleration algorithms
and provides a type-safe dispatch system to replace the previous symbol-based approach.
"""

# Abstract type hierarchy
abstract type FixedPointAlgorithm end

# Simple iteration (no acceleration)
"""
    Simple()

Simple fixed point iteration: x_{n+1} = f(x_n)

This is the most basic fixed point iteration method with no acceleration.
It's guaranteed to converge for contractive functions but may be slow.
"""
struct Simple <: FixedPointAlgorithm end

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
        condition_threshold >= 1.0 || throw(ArgumentError("condition_threshold must be at least 1.0"))
        new(maxM, condition_threshold)
    end
end

"""
    Aitken()

Aitken's Δ² acceleration method.

This scalar acceleration method (applied elementwise to vectors) assumes the sequence
converges at a constant rate and uses three consecutive iterates to estimate the limit.

# References
Aitken, A.C. (1926). "On Bernoulli's Numerical Solution of Algebraic Equations."
Proceedings of the Royal Society of Edinburgh 46: 289-305.
"""
struct Aitken <: FixedPointAlgorithm end

"""
    Newton()

Newton-based acceleration using approximate derivatives.

This method uses finite differences to approximate the derivative and applies
Newton's method to solve f(x) - x = 0.
"""
struct Newton <: FixedPointAlgorithm end

# Polynomial extrapolation methods
"""
    MPE(; extrapolation_period=7)

Minimal Polynomial Extrapolation.

MPE takes a linear combination of previous sequential iterates to accelerate convergence.
Extrapolation is performed every `extrapolation_period` iterations.

# Parameters
- `extrapolation_period::Int`: Number of simple iterations before extrapolation (default: 7)

# References
Cabay, S., and L.W. Jackson (1976). "A Polynomial Extrapolation Method for Finding
Limits and Antilimits of Vector Sequences." SIAM Journal of Numerical Analysis 13(5): 734-52.
"""
struct MPE <: FixedPointAlgorithm
    extrapolation_period::Int

    function MPE(; extrapolation_period::Int=7)
        extrapolation_period >= 1 || throw(ArgumentError("extrapolation_period must be at least 1"))
        new(extrapolation_period)
    end
end

"""
    RRE(; extrapolation_period=7)

Reduced Rank Extrapolation.

RRE is similar to MPE but uses a different linear combination approach.
Extrapolation is performed every `extrapolation_period` iterations.

# Parameters
- `extrapolation_period::Int`: Number of simple iterations before extrapolation (default: 7)

# References
Eddy, R.P. (1979). "Extrapolating to the Limit of a Vector Sequence."
Information Linkage Between Applied Mathematics and Industry: 387-96.
"""
struct RRE <: FixedPointAlgorithm
    extrapolation_period::Int

    function RRE(; extrapolation_period::Int=7)
        extrapolation_period >= 1 || throw(ArgumentError("extrapolation_period must be at least 1"))
        new(extrapolation_period)
    end
end

# Epsilon algorithms
"""
    VEA(; extrapolation_period=7)

Vector Epsilon Algorithm.

VEA uses the epsilon algorithm with Moore-Penrose pseudoinverse for vector problems.
Extrapolation is performed every `extrapolation_period` iterations.

# Parameters
- `extrapolation_period::Int`: Number of simple iterations before extrapolation (default: 7)
                              Should be even for optimal performance.

# References
Wynn, P. (1956). "On a Device for Computing the e_m(S_n) Transformation."
Mathematical Tables and Other Aids to Computation 10(54): 91-96.
"""
struct VEA <: FixedPointAlgorithm
    extrapolation_period::Int

    function VEA(; extrapolation_period::Int=7)
        extrapolation_period >= 1 || throw(ArgumentError("extrapolation_period must be at least 1"))
        new(extrapolation_period)
    end
end

"""
    SEA(; extrapolation_period=7)

Scalar Epsilon Algorithm.

SEA applies the epsilon algorithm elementwise to vector problems.
Extrapolation is performed every `extrapolation_period` iterations.

# Parameters
- `extrapolation_period::Int`: Number of simple iterations before extrapolation (default: 7)
                              Should be even for optimal performance.

# References
Wynn, P. (1956). "On a Device for Computing the e_m(S_n) Transformation."
Mathematical Tables and Other Aids to Computation 10(54): 91-96.
"""
struct SEA <: FixedPointAlgorithm
    extrapolation_period::Int

    function SEA(; extrapolation_period::Int=7)
        extrapolation_period >= 1 || throw(ArgumentError("extrapolation_period must be at least 1"))
        new(extrapolation_period)
    end
end

# Utility functions for the algorithm types

"""
    algorithm_name(alg::FixedPointAlgorithm)

Get the name of the algorithm as a string for display purposes.
"""
algorithm_name(::Simple) = "Simple"
algorithm_name(::Anderson) = "Anderson"
algorithm_name(::Aitken) = "Aitken"
algorithm_name(::Newton) = "Newton"
algorithm_name(::MPE) = "MPE"
algorithm_name(::RRE) = "RRE"
algorithm_name(::VEA) = "VEA"
algorithm_name(::SEA) = "SEA"

"""
    needs_extrapolation_period(alg::FixedPointAlgorithm)

Check if the algorithm uses an extrapolation period.
"""
needs_extrapolation_period(::Union{Simple,Anderson,Aitken,Newton}) = false
needs_extrapolation_period(::Union{MPE,RRE,VEA,SEA}) = true

"""
    get_extrapolation_period(alg::FixedPointAlgorithm)

Get the extrapolation period for algorithms that use it.
"""
get_extrapolation_period(alg::Union{MPE,RRE,VEA,SEA}) = alg.extrapolation_period
get_extrapolation_period(::Union{Simple,Anderson,Aitken,Newton}) = 1  # Default for non-extrapolation methods

"""
    is_polynomial_method(alg::FixedPointAlgorithm)

Check if the algorithm uses polynomial extrapolation.
"""
is_polynomial_method(::Union{MPE,RRE}) = true
is_polynomial_method(::FixedPointAlgorithm) = false

"""
    is_epsilon_method(alg::FixedPointAlgorithm)

Check if the algorithm uses epsilon extrapolation.
"""
is_epsilon_method(::Union{VEA,SEA}) = true
is_epsilon_method(::FixedPointAlgorithm) = false

# Legacy symbol support (for internal use with extrapolation methods)
# These are needed for PolynomialExtrapolation and EpsilonExtrapolation functions
algorithm_to_symbol(::Simple) = :Simple
algorithm_to_symbol(::Anderson) = :Anderson
algorithm_to_symbol(::Aitken) = :Aitken
algorithm_to_symbol(::Newton) = :Newton
algorithm_to_symbol(::MPE) = :MPE
algorithm_to_symbol(::RRE) = :RRE
algorithm_to_symbol(::VEA) = :VEA
algorithm_to_symbol(::SEA) = :SEA

# Algorithm implementations are now in algorithm_implementations.jl to avoid method ambiguities
