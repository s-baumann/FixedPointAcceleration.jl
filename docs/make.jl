using Distributions, Documenter, FixedPointAcceleration

makedocs(
    format = :html,
    sitename = "FixedPointAcceleration",
    modules = [FixedPointAcceleration]
)

deploydocs(
    repo   = "github.com/s-baumann/FixedPointAcceleration.jl.git",
    julia  = "1.0",
    target = "build",
    deps   = nothing,
    make   = nothing
)
