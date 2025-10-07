module test_CommonIncrements

using Test
@testset "Test Common Increments" begin
    using FixedPointAcceleration
    function funcfunc(x::Array{Float64,1})
        # The first coordinate convergences to 4.0 by 1 unit per iterate.
        output = Array{Float64,1}(x)
        if abs(x[1] - 4.0) <= 1.0
            output[1] = 4.0
        elseif x[1] > 4.0
            output[1] = x[1] - 1.0
        else
            output[1] = x[1] + 1.0
        end
        # The second does aitken convergence to 2.3
        output[2] = x[2] + (2.3-x[2])/2.0
        return output
    end
    Inputs = [19.0, 10.0]
    fp_anderson = fixed_point(
        funcfunc,
        Inputs,
        Anderson(),
        FixedPointOptions(quiet_errors=true, replace_invalids=:ReplaceElements),
    )
    @test fp_anderson.termination_condition == :ReachedConvergenceThreshold
    fp_aitken = fixed_point(
        funcfunc,
        Inputs,
        Aitken(),
        FixedPointOptions(quiet_errors=true, replace_invalids=:ReplaceElements),
    )
    @test fp_aitken.termination_condition == :ReachedConvergenceThreshold
    fp_newton = fixed_point(
        funcfunc,
        Inputs,
        Newton(),
        FixedPointOptions(quiet_errors=true, replace_invalids=:ReplaceElements),
    )
    @test fp_newton.termination_condition == :ReachedConvergenceThreshold
    fp_simple = fixed_point(
        funcfunc,
        Inputs,
        Simple(),
        FixedPointOptions(
            quiet_errors=true, replace_invalids=:ReplaceElements, reporting_sig_figs=10
        ),
    )
    @test fp_simple.termination_condition == :ReachedConvergenceThreshold
    fp_SEA = fixed_point(
        funcfunc,
        Inputs,
        SEA(),
        FixedPointOptions(quiet_errors=true, replace_invalids=:ReplaceElements),
    )
    @test fp_SEA.termination_condition == :ReachedConvergenceThreshold
    fp_VEA = fixed_point(
        funcfunc,
        Inputs,
        VEA(),
        FixedPointOptions(quiet_errors=true, replace_invalids=:ReplaceElements),
    )
    @test fp_VEA.termination_condition == :ReachedConvergenceThreshold
    fp_RRE = fixed_point(
        funcfunc,
        Inputs,
        RRE(),
        FixedPointOptions(quiet_errors=true, replace_invalids=:ReplaceElements),
    )
    @test fp_RRE.termination_condition == :ReachedMaxIter # This one fails because it keeps proposing bad ideas.
    fp_MPE = fixed_point(
        funcfunc,
        Inputs,
        MPE(),
        FixedPointOptions(quiet_errors=true, replace_invalids=:ReplaceElements),
    )
    @test fp_MPE.termination_condition == :ReachedConvergenceThreshold
end

end
