---
title: 'FixedPointAcceleration: A Julia package for accelerating convergence of fixed point problems'
tags:
  - Julia
  - Fixed Point Acceleration
  - dynamics
authors:
  - name: Stuart Baumann
    orcid: 0000-0002-9657-0969
    affiliation: Credit Suisse Securities (Europe)
  - name: Margaryta Klymak
    orcid: 0000-0003-4376-883X
    affiliation: University of Oxford
affiliations:
 - name: Credit Suisse Securities (Europe)
   index: 1
 - name: University of Oxford
   index: 2
date: 12 December 2020
bibliography: paper.bib

---

# Summary

A fixedpoint problem is one where we have a function taking and returning a function, $f: X^d -> X^d$ for which we seek a value $x$ such that $f(x) = x$.
The simplest way to solve these problems is to choose an initial guess, $g$ and then estimate the sequence f(g), f(f(g)), f(f(f(g))), and so on.
In many cases this sequence will converge to the fixed point.
While this method is robust it is also typically slow. 
FixedPointAcceleration is a Julia package which implements a number of extrapolation 
algorithms to accelerate the finding of a fixedpoint.
These algorithms include Newton acceleration, Aitken acceleration and Anderson acceleration [@Anderson1965;@Walker2010;@WalkerNi2011] as well as epsilon extrapolation methods [@Wynn1962] and minimal polynomial methods [@CabayJackson1976;@SmithFordSidi1987].

There are a large number of research applications for fixed point acceleration.
In statistics they can be used to speed up the estimation of Expectation Maximisation (EM) algorithms.
There are also applications in consumption smoothing problems in economics.
The package documentation presents examples of how FixedPointAcceleration can be used to solve problems in these two areas as well as other examples in economics, statistics and finance.

Prior to the release of FixedPointAcceleration there were no Julia packages offering fixed point acceleration.
There was however an earlier package in R (``FixedPoint'') by the same authors that implements almost identical functionality as FixedPointAcceleration.
This R package is documented in @BaumannKlymak2019 .
While Julia is designed to facilitate calling packages from other languages, there were a few reasons why the reimplementation of these algorithms was preferable.
The first is that as a compiled language Julia will run faster.
In practical terms this means that the use of a fixed point acceleration algorithm will be useful for finding fixedpoints of cheap functions as well as expensive ones.
The second is that Julia allows automatic differentiation but only if all code being executed is Julia code.
The implementation of fixed point algorithms in Julia allows it to be used in code without preventing the passage of gradiants.

# References