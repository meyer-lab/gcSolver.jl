__precompile__()
module gcSolver

using OrdinaryDiffEq
import LinearAlgebra: diag
import ForwardDiff
using Optim
using Statistics

include("reaction.jl")

const solTol = 1.0e-9

function domainDef(u, p, t)
    return any(x -> x < -solTol, u)
end


" Check that we haven't been provided anything ridiculous by the user. "
function checkInputs(tps::Vector, params::Vector)
    @assert all(tps .>= 0.0)
    @assert length(params) == Nparams
    @assert all(params .>= 0.0)
    @assert params[25] < 1.0
end


function runCkineSetup(tps::Vector{Float64}, params::Vector)
    checkInputs(tps, params)
    u0 = solveAutocrine(params)

    return ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), params)
end


" Actually run the gc ODE model. "
function runCkine(tps::Vector{Float64}, params::Vector)::Matrix
    prob = runCkineSetup(tps, params)

    alg = AutoTsit5(Rodas5(autodiff = false))
    sol = solve(prob, alg; saveat = tps, reltol = solTol, isoutofdomain = domainDef).u

    if length(tps) > 1
        sol = vcat(transpose.(sol)...)
    else
        sol = reshape(sol[1], (1, Nspecies))
    end

    return sol
end


" Converts the ODE solution to a predicted amount of pSTAT. "
function runCkinePSTAT(tps::Vector, params::Vector)::Vector
    sol = runCkine(tps, params)

    # Summation of active species
    pSTAT = sol[:, 43] + 2 * (sol[:, 44] + sol[:, 45])

    @assert length(pSTAT) == length(tps)
    return vec(pSTAT)
end


" Calculate the Jacobian of the model and perform variance propagation. "
function runCkineVarProp(tps::Vector, params::Vector, sigma)::Vector
    checkInputs(tps, params)

    # Sigma is the covariance matrix of the input parameters
    function jacF(x)
        pp = vcat(params[1:27], x, params[33:end])
        return runCkinePSTAT(tps, pp)
    end

    jac = zeros(5, length(tps))
    ForwardDiff.jacobian!(jac, jacF, params[28:32])

    # Just return the diagonal for the marginal variance
    return diag(transpose(jac) * sigma * jac)
end


export runCkine, runCkineVarProp

end # module
