using Documenter, FixedPointAcceleration

makedocs(
    format = Documenter.HTML(),
    sitename = "FixedPointAcceleration",
    modules = [FixedPointAcceleration]
)

deploydocs(
    repo   = "github.com/s-baumann/FixedPointAcceleration.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
