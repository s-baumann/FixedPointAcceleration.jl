module test_quality_assurance

using Test
@testset "Aqua" begin
    using FixedPointAcceleration
    using Aqua: Aqua
    Aqua.test_all(FixedPointAcceleration)
end

end
