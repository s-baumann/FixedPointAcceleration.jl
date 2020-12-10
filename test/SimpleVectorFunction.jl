func(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
Inputs = [0.3,900.0]
fp_simple   = fixed_point(func, Inputs; Algorithm = :Simple)
fp_simple.Convergence_ < 1e-10
fp_anderson = fixed_point(func, Inputs; Algorithm = :Anderson)
fp_anderson.Convergence_ < 1e-10
fp_aitken   = fixed_point(func, Inputs; Algorithm = :Aitken)
fp_aitken.Convergence_ < 1e-10
fp_newton   = fixed_point(func, Inputs; Algorithm = :Newton)
fp_newton.Convergence_ < 1e-10
fp_VEA      = fixed_point(func, Inputs; Algorithm = :VEA)
fp_VEA.Convergence_ < 1e-10
fp_SEA      = fixed_point(func, Inputs; Algorithm = :SEA)
fp_SEA.Convergence_ < 1e-10
fp_MPE      = fixed_point(func, Inputs; Algorithm = :MPE)
fp_MPE.Convergence_ < 1e-10
fp_RRE      = fixed_point(func, Inputs; Algorithm = :RRE)
fp_RRE.Convergence_ < 1e-10

# Now trying a function with a typeswitch from Int to Float
Inputs = [1,4]
fp_simple   = fixed_point(func, Inputs; Algorithm = :Simple)
fp_simple.Convergence_ < 1e-10
fp_anderson = fixed_point(func, Inputs; Algorithm = :Anderson)
fp_anderson.Convergence_ < 1e-10

# And finally printing a couple of them to test printing.
fp_anderson2 = fixed_point(func, Inputs; Algorithm = :Anderson, PrintReports = true)
fp_anderson2.Convergence_ < 1e-10
fp_simple2   = fixed_point(func, Inputs; Algorithm = :Simple, PrintReports = true)
fp_simple2.Convergence_ < 1e-10

# Testing function that returns a missing
f(x) = [missing, missing]
fp_simple3   = fixed_point(f, Inputs; Algorithm = :Simple, PrintReports = true)
fp_simple3.TerminationCondition_ == :InvalidInputOrOutputOfIteration
f(x) = [1, missing]
fp_simple4   = fixed_point(f, Inputs; Algorithm = :Simple, PrintReports = true)
fp_simple4.TerminationCondition_ == :InvalidInputOrOutputOfIteration

# Testing the outputting of a tuple
f(x) = (a = [1.0, 1.0], b = :goodProgressInFunction)
@test_throws ErrorException("The Fixedpoint function can only return a vector or a tuple  of which the first entry is the vector for which a fixedpoint is sought and the second is a namedtuple (the contents of which are output for the user but are not used in fixed point acceleration).") fixed_point(f, Inputs; Algorithm = :Simple, PrintReports = true)
# And doing side effects properly.
f(x) = ([1.0, 1.0], (b = :goodProgressInFunction,))
fp_simple5 = fixed_point(f, Inputs; Algorithm = :Simple, PrintReports = true)
fp_simple5.Convergence_ < 1e-10
