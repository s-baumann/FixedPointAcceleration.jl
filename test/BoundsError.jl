using FixedPointAccelerationNext, Test

@testset "Test bounds" begin
    # Testing Error Evaluating Function
    simple_vector_function(x) = [0.5*sqrt(x[1] + x[2]), 1.5*x[1] + 0.5*x[2]]
    Inputs = [0.3, 900]
    opts_quiet = FixedPointOptions(quiet_errors=true)
    fp = fixed_point(simple_vector_function, Inputs, Anderson(), opts_quiet)
    # Inspecting this fp reveals an error after the 3rd iteration because
    # Anderson tries to use a negative value for both x entries which results in
    # the square root of a negative number. We can switch to simple
    # iterations for a while to fix this.
    @test fp.termination_condition == :InvalidInputOrOutputOfIteration
    @test fp.failed_evaluation.error == :ErrorExecutingFunction
    opts_7 = FixedPointOptions(max_iterations=7)
    fp = fixed_point(simple_vector_function, fp, Simple(), opts_7)
    @test fp.termination_condition == :ReachedMaxIter
    fp = fixed_point(simple_vector_function, fp, Anderson())
    @test fp.termination_condition == :ReachedConvergenceThreshold

    # Testing Input of NaN
    fp = fixed_point(simple_vector_function, [NaN, 900], Simple(), opts_quiet)
    @test fp.failed_evaluation.error == :InputNAsDetected
    # Testing Input of Inf
    fp = fixed_point(simple_vector_function, [-Inf, 900], Simple(), opts_quiet)
    @test fp.failed_evaluation.error == :InputInfsDetected

    # Testing Output of Nan
    function funcfunc1(x::Array{Float64,1})
        if abs(x[1] - 4.0) < 1e-12
            return Array{Float64,1}([NaN, 4.0])
        end
        return sqrt.(x)
    end
    Inputs = [4.0, 1.0]
    fp = fixed_point(funcfunc1, Inputs, Anderson(), opts_quiet)
    @test fp.failed_evaluation.error == :OutputNAsDetected
    # Testing Output of Missing
    function funcfunc2(x::Array{Float64,1})
        if abs(x[1] - 4.0) < 1e-12
            return [missing, 4.0]
        end
        return sqrt.(x)
    end
    Inputs = [4.0, 1.0]
    fp = fixed_point(funcfunc2, Inputs, Anderson(), opts_quiet)
    @test fp.failed_evaluation.error == :OutputMissingsDetected
    # Testing Output of Inf
    function funcfunc3(x::Array{Float64,1})
        if abs(x[1] - 4.0) < 1e-12
            return Array{Float64,1}([Inf, 4.0])
        end
        return sqrt.(x)
    end
    Inputs = [4.0, 1.0]
    fp = fixed_point(funcfunc3, Inputs, Anderson(), opts_quiet)
    @test fp.failed_evaluation.error == :OutputInfsDetected
    # Testing Output of wrong size
    function funcfunc4(x::Array{Float64,1})
        if abs(x[1] - 4.0) < 1e-12
            return Array{Float64,1}([5.0, 4.0, 4.0])
        end
        return sqrt.(x)
    end
    Inputs = [4.0, 1.0]
    fp = fixed_point(funcfunc4, Inputs, Anderson(), opts_quiet)
    @test fp.failed_evaluation.error == :LengthOfOutputNotSameAsInput

    # Testing Output of wrong type
    function funcfunc5(x::Array{Float64,1})
        return Array{Int,1}([5, 4])
    end
    Inputs = [4.0, 1.0]
    fp = fixed_point(funcfunc5, Inputs, Anderson(), opts_quiet)
    @test fp.failed_evaluation.error == :FunctionIsNotTypeStable
end
