using Test, FixedPointAccelerationNext

@testset "Aqua" begin
    using Aqua: Aqua
    Aqua.test_all(FixedPointAccelerationNext)
end

@testset "JET" begin
    using JET
    JET.test_package(FixedPointAccelerationNext; target_defined_modules=true)
end
