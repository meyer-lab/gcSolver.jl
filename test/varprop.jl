
@testset "Reasonable return from varprop function." begin
    sigma = Matrix{Int}(I, 5, 5)

    retval = runCkineVarProp(tps, rxntfR, sigma)

    @test length(retval) == length(tps)
    @test all(retval .>= 0.0)
end


@testset "Test cost gradient function." begin
    g = ForwardDiff.gradient(x -> runCkineCost(tps, x, tps), rxntfR)

    @time ForwardDiff.gradient(x -> runCkineCost(tps, x, tps), rxntfR)

    @test length(g) == length(rxntfR)
    @test all(isfinite.(g))
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
