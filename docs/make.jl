using Documenter, FixedPointAcceleration

makedocs(
    format = Documenter.HTML(),
    sitename = "FixedPointAcceleration",
    modules = [FixedPointAcceleration],
    pages = ["Index" => "index.md",
             "Introduction" => "1_FixedPoints.md",
             "Algorithms" => "2_Algorithms.md",
             "Advice in using package" => "3_UsingAdvice.md",
             "Applications" => "4_Applications.md",
             "Termination Conditions" => "5_TerminationConditions.md",
             "API" => "api.md"
             "References" => "99_refs.md"]
)

deploydocs(
    repo   = "github.com/s-baumann/FixedPointAcceleration.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
