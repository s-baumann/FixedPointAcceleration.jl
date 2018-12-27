func(x) = cos.(x)
Inputs = 1.1
fp_simple   = fixed_point(func, Inputs; Algorithm = Simple)
fp_simple.Convergence_ < 1e-10
fp_anderson = fixed_point(func, Inputs; Algorithm = Anderson)
fp_anderson.Convergence_ < 1e-10
fp_aitken   = fixed_point(func, Inputs; Algorithm = Aitken)
fp_aitken.Convergence_ < 1e-10
fp_newton   = fixed_point(func, Inputs; Algorithm = Newton)
fp_newton.Convergence_ < 1e-10
fp_VEA      = fixed_point(func, Inputs; Algorithm = VEA)
fp_VEA.Convergence_ < 1e-10
fp_SEA      = fixed_point(func, Inputs; Algorithm = SEA)
fp_SEA.Convergence_ < 1e-10
fp_MPE      = fixed_point(func, Inputs; Algorithm = MPE)
fp_MPE.Convergence_ < 1e-10
fp_RRE      = fixed_point(func, Inputs; Algorithm = RRE)
fp_RRE.Convergence_ < 1e-10
