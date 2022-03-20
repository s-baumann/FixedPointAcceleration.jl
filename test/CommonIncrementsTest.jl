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
    Inputs = [19.0,10.0]
    fp_anderson = fixed_point(funcfunc, Inputs; Algorithm = :Anderson, PrintReports = false, ReplaceInvalids = :ReplaceElements)
    @test fp_anderson.TerminationCondition_ == :ReachedConvergenceThreshold
    fp_aitken = fixed_point(funcfunc, Inputs; Algorithm = :Aitken, PrintReports = false, ReplaceInvalids = :ReplaceElements)
    @test fp_aitken.TerminationCondition_ == :ReachedConvergenceThreshold
    fp_newton = fixed_point(funcfunc, Inputs; Algorithm = :Newton, PrintReports = false, ReplaceInvalids = :ReplaceElements)
    @test fp_newton.TerminationCondition_ == :ReachedConvergenceThreshold
    fp_simple = fixed_point(funcfunc, Inputs; Algorithm = :Simple, PrintReports = false, ReplaceInvalids = :ReplaceElements, ReportingSigFig =  Integer(10))
    @test fp_simple.TerminationCondition_ == :ReachedConvergenceThreshold
    fp_SEA = fixed_point(funcfunc, Inputs; Algorithm = :SEA, PrintReports = false, ReplaceInvalids = :ReplaceElements)
    @test fp_SEA.TerminationCondition_ == :ReachedConvergenceThreshold
    fp_VEA = fixed_point(funcfunc, Inputs; Algorithm = :VEA, PrintReports = false, ReplaceInvalids = :ReplaceElements)
    @test fp_VEA.TerminationCondition_ == :ReachedConvergenceThreshold
    fp_RRE = fixed_point(funcfunc, Inputs; Algorithm = :RRE, PrintReports = false, ReplaceInvalids = :ReplaceElements)
    @test fp_RRE.TerminationCondition_ == :ReachedMaxIter # This one fails because it keeps proposing bad ideas.
    fp_MPE = fixed_point(funcfunc, Inputs; Algorithm = :MPE, PrintReports = false, ReplaceInvalids = :ReplaceElements)
    @test fp_MPE.TerminationCondition_ == :ReachedConvergenceThreshold
end
