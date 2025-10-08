using Test, FixedPointAccelerationNext

@testset "Aqua" begin
    using Aqua: Aqua
    Aqua.test_all(FixedPointAccelerationNext)
end

@testset "JET" begin
    using JET
    JET.test_package(FixedPointAccelerationNext; target_defined_modules=true)
end

@testset "ExplicitImports" begin
    using ExplicitImports

    @test check_no_implicit_imports(FixedPointAccelerationNext) == nothing
    @test check_all_explicit_imports_via_owners(FixedPointAccelerationNext) == nothing
    @test check_all_explicit_imports_are_public(FixedPointAccelerationNext) == nothing
    @test check_no_stale_explicit_imports(FixedPointAccelerationNext) == nothing
    @test check_all_qualified_accesses_via_owners(FixedPointAccelerationNext) == nothing
    @test check_all_qualified_accesses_are_public(FixedPointAccelerationNext) == nothing
    @test check_no_self_qualified_accesses(FixedPointAccelerationNext) == nothing
end

@testset "CheckConcreteStructs" begin
    using CheckConcreteStructs
    all_concrete(FixedPointAccelerationNext.FixedPointConfig)
    all_concrete(FixedPointAccelerationNext.FixedPointSolution)

    all_concrete(FixedPointAccelerationNext.Simple)
    all_concrete(FixedPointAccelerationNext.Anderson)
    all_concrete(FixedPointAccelerationNext.Aitken)
    all_concrete(FixedPointAccelerationNext.MPE)
    all_concrete(FixedPointAccelerationNext.RRE)
    all_concrete(FixedPointAccelerationNext.VEA)
    all_concrete(FixedPointAccelerationNext.SEA)

    all_concrete(FixedPointAccelerationNext.ProgressTracker)
    all_concrete(FixedPointAccelerationNext.IterationState)
    all_concrete(FixedPointAccelerationNext.Workspace)
    all_concrete(FixedPointAccelerationNext.IterationCallbacks)
end
