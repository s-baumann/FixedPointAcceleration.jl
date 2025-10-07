"""
Base types and utilities for fixed point acceleration algorithms.

This module defines the abstract type hierarchy and common utilities
used by all fixed point acceleration algorithms.
"""

# Abstract type hierarchy
abstract type FixedPointAlgorithm end

"""
    algorithm_name(alg::FixedPointAlgorithm)

Get the name of the algorithm as a string for display purposes.
"""
function algorithm_name end

"""
    needs_extrapolation_period(alg::FixedPointAlgorithm)

Check if the algorithm uses an extrapolation period.
"""
function needs_extrapolation_period end

"""
    get_extrapolation_period(alg::FixedPointAlgorithm)

Get the extrapolation period for algorithms that use it.
"""
function get_extrapolation_period end

"""
    is_polynomial_method(alg::FixedPointAlgorithm)

Check if the algorithm uses polynomial extrapolation.
"""
function is_polynomial_method end

"""
    is_epsilon_method(alg::FixedPointAlgorithm)

Check if the algorithm uses epsilon extrapolation.
"""
function is_epsilon_method end

# Default implementations
needs_extrapolation_period(::FixedPointAlgorithm) = false
get_extrapolation_period(::FixedPointAlgorithm) = 1
is_polynomial_method(::FixedPointAlgorithm) = false
is_epsilon_method(::FixedPointAlgorithm) = false
