using Test
@testset "Aqua" begin
    using FixedPointAcceleration
    using Aqua: Aqua
    Aqua.test_all(FixedPointAcceleration)
end
