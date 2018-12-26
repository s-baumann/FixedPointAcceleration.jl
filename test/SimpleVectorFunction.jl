func(x) = [0.5*sqrt(abs(x[1] + x[2])), 1.5*x[1] + 0.5*x[2]]
Inputs = [1.1,2.2]
fp_simple   = fixed_point(func, Inputs; Algorithm = Simple)
fp_simple.Convergence_ < 1e-10
fp_anderson = fixed_point(func, Inputs; Algorithm = Anderson)
fp_aitken   = fixed_point(func, Inputs; Algorithm = Aitken)
fp_newton   = fixed_point(func, Inputs; Algorithm = Newton)
fp_VEA      = fixed_point(func, Inputs; Algorithm = VEA)
fp_SEA      = fixed_point(func, Inputs; Algorithm = SEA)
fp_MPE      = fixed_point(func, Inputs; Algorithm = MPE)
fp_RRE      = fixed_point(func, Inputs; Algorithm = RRE)
fp_RRE.Convergence_ < 1e-10
