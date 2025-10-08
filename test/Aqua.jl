using Test, FixedPointAcceleration

@testset "Aqua" begin
    using Aqua: Aqua
    Aqua.test_all(FixedPointAcceleration)
end
