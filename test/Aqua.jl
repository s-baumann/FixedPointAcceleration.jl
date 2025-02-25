using Test
@testset "Aqua" begin
    using FixedPointAcceleration
    import Aqua
    Aqua.test_all(FixedPointAcceleration)
end
