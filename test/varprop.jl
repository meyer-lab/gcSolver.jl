
@testset "Reasonable return from varprop function." begin
    sigma = Matrix{Int}(I, 3, 3)
    recAbund = ones(5)

    retval = runCkineVarProp(tps, rxntfR, sigma, recAbund)

    @test length(retval) == length(tps)
    @test all(retval .>= 0.0)
end


@testset "Run fitting." begin
    p = gcSolver.getUnkVec()
    outtG = similar(p)

    outt = gcSolver.resids(p)
    ForwardDiff.gradient!(outtG, gcSolver.resids, p)

    @time ForwardDiff.gradient!(outtG, gcSolver.resids, p)

    @test isfinite(outt)
    @test all(isfinite.(outtG))
end


@testset "Test that residuals can be calculated using Farhat Fit." begin
    farhatVec = gcSolver.getFarhatVec()
    resids = gcSolver.resids(farhatVec)
    @test isfinite(resids)
end
