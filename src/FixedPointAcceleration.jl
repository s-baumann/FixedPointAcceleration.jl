module FixedPointAcceleration
using Dates: now

using LinearAlgebra: qr, norm, pinv, qr!, QRPivoted, norm, svd, Diagonal

include("config.jl")
include("state.jl")
include("methods.jl")
include("workspace.jl")
include("anderson_impl.jl")
include("polynomial_impl.jl")
include("epsilon_impl.jl")
include("solve.jl")
include("step.jl")

export solve, solve!, FixedPointConfig, FixedPointSolution
export Simple, Anderson, Aitken, MPE, RRE, VEA, SEA

include("old/FixedPointAcceleration.jl")
export OldImplementation

end
