
function gradFunc(x)
	return runCkine(tps, x)
end


@testset "Profile forward sensitivities." begin
	jac = zeros(gcSolver.Nspecies*length(tps), gcSolver.Nparams)

    ForwardDiff.jacobian!(jac, gradFunc, rxntfR)

    @time ForwardDiff.jacobian!(jac, gradFunc, rxntfR)
end
