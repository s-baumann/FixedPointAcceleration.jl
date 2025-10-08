using Test, FixedPointAccelerationNext

@testset "Aqua" begin
    using Aqua: Aqua
    Aqua.test_all(FixedPointAccelerationNext)
end
