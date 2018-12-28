func(x) = cos.(x)
Inputs = 1.1
# We should not have any slowdown from starting in simple for two. Then VEA for 3. Then Simple. Then VEA until the end.
# The reason is that VEA does nothing of interest until 7 by default anyway.
fp_chain      = fixed_point(func, Inputs; Algorithm = Simple, MaxIter = 2)
fp_chain      = fixed_point(func, fp_chain; Algorithm = VEA, MaxIter = 3)
fp_chain      = fixed_point(func, fp_chain; Algorithm = Simple, MaxIter = 1)
fp_chain      = fixed_point(func, fp_chain; Algorithm = VEA, MaxIter = 100)

fp_nochain = fixed_point(func, Inputs; Algorithm = VEA, MaxIter = 100)
fp_chain.Iterations_ == fp_nochain.Iterations_
all(abs.(fp_nochain.Inputs_ .- fp_chain.Inputs_) .< 1e-14)
