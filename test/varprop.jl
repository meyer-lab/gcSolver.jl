
function gradFunc(x)
	return runCkine(tps, x)
end


@testset "Profile forward sensitivities." begin
	jac = zeros(gcSolver.Nspecies*length(tps), gcSolver.Nparams)

    ForwardDiff.jacobian!(jac, gradFunc, rxntfR)

    @time ForwardDiff.jacobian!(jac, gradFunc, rxntfR)
end


@testset "Reasonable return from varprop function." begin
    sigma = Matrix{Int}(I, length(rxntfR), length(rxntfR))

    retval = gcSolver.runCkineVarPorp(tps, rxntfR, sigma)

    @time gcSolver.runCkineVarPorp(tps, rxntfR, sigma)

    @test size(retval) == (length(tps), length(tps))
end

