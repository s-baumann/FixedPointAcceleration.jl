module test_replace_invalids

using Test

using FixedPointAcceleration

@testset "Replace Invalids Testing" begin
    # Testing algos generating an invalid input.
    function funcfunc(x::Array{Float64,1})
        # The first coordinate convergence to 4.0 by 1 unit per iterate.
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
    fp = fixed_point(funcfunc, Inputs, Aitken())
    @test fp.failed_evaluation.error == :InputInfsDetected
    @test fp.failed_evaluation.input[1] == -Inf
    @test !isinf(fp.failed_evaluation.input[2])
    # Now fixing with replace element
    opts_replace = FixedPointOptions(replace_invalids=:ReplaceElements)
    fp = fixed_point(funcfunc, Inputs, Aitken(), opts_replace)
    @test fp.termination_condition == :ReachedConvergenceThreshold
    # Now fixing with replace element
    opts_replace_vector = FixedPointOptions(replace_invalids=:ReplaceVector)
    fp2 = fixed_point(funcfunc, Inputs, Aitken(), opts_replace_vector)
    @test fp2.termination_condition == :ReachedConvergenceThreshold
    @test fp.iterations != fp2.iterations
end

end
