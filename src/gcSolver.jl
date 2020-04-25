module gcSolver

using OrdinaryDiffEq
import LinearAlgebra: diag, norm
import ForwardDiff
using Optim
using Statistics
import ModelingToolkit

include("reaction.jl")

const solTol = 1.0e-12

function domainDef(u, p, t)
    return any(x -> x < -solTol, u)
end


" Check that we haven't been provided anything ridiculous by the user. "
function checkInputs(tps::Vector{Float64}, params::Vector)
    @assert all(tps .>= 0.0)
    @assert length(params) == Nparams
    @assert all(params .>= 0.0)
    @assert params[22] < 1.0
end


" This recompiles the ODE function into a symbolic Jacobian. "
function modelCompile()
    u0 = ones(Nspecies)
    params = ones(Nparams) * 0.1

    prob = ODEProblem(fullDeriv, u0, (0.0, 1.0), params)
    deMT, varsMT, paramsMT = ModelingToolkit.modelingtoolkitize(prob)

    f_iip = eval(ModelingToolkit.generate_function(deMT, varsMT, paramsMT)[2])
    tgrad_iip = eval(ModelingToolkit.generate_tgrad(deMT)[2])
    jac = eval(ModelingToolkit.generate_jacobian(deMT)[2])

    return ODEFunction(f_iip; tgrad = tgrad_iip, jac = jac)
end


const modelFunc = modelCompile()


function runCkineSetup(tps::Vector{Float64}, params::Vector)
    checkInputs(tps, params)
    u0 = solveAutocrine(params)

    # Remove receptor expression if we're modeling no-trafficking
    if params[20] == 0.0
        params[25:29] .= 0.0  # set receptor expression to 0.0
    end

    return ODEProblem(modelFunc, u0, (0.0, maximum(tps)), params)
end


" Actually run the gc ODE model. "
function runCkine(tps::Vector{Float64}, params::Vector)::Matrix
    prob = runCkineSetup(tps, params)

    sol = solve(prob, AutoTsit5(Rodas5(); nonstifftol = 10 // 10); saveat = tps, reltol = solTol, isoutofdomain = domainDef).u

    if length(tps) > 1
        sol = vcat(transpose.(sol)...)
    else
        sol = reshape(sol[1], (1, Nspecies))
    end

    return sol
end


" Converts the ODE solution to a predicted amount of pSTAT. "
function runCkinePSTAT(tps::Vector{Float64}, params::Vector)::Vector
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

include("fit.jl")

export runCkine, runCkineVarProp

end # module
