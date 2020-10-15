module gcSolver

using OrdinaryDiffEq
import LinearAlgebra: diag, norm, dot
import ForwardDiff
import ForwardDiff: value, Dual, partials, jacobian!, jacobian
using Optim
using Statistics
import ModelingToolkit: modelingtoolkitize
import CSV
using DataFrames
import StatsBase: indicatormat
using StatsFuns

include("reaction.jl")
include("dataImport.jl")


" Check that we haven't been provided anything ridiculous by the user. "
function checkInputs(tps::Vector{Float64}, params::Vector)
    @assert all(tps .>= 0.0)
    @assert params[22] < 1.0
    @assert length(params) == Nparams
    @assert all(params .>= 0.0)
end


" This recompiles the ODE function into a symbolic Jacobian. "
function __init__()
    u0 = ones(Nspecies)
    params = 0.1ones(Nparams)

    prob = ODEProblem(fullDeriv, u0, (0.0, 1.0), params)
    deMT = modelingtoolkitize(prob)

    global modelFunc = ODEFunction(deMT; jac = true)
end


function runCkineSetup(tps::Vector{Float64}, p::Vector{T}) where {T}
    checkInputs(tps, p)
    u0 = solveAutocrine(p)

    return ODEProblem(modelFunc, u0, convert(T, maximum(tps)), p)
end


" Actually run the gc ODE model. "
function runCkine(tps::Vector{Float64}, params; pSTAT5 = false)
    if params isa Vector
        prob = runCkineSetup(tps, params)
    else
        prob = params
    end

    if pSTAT5
        sidx = pSTATidx
    else
        sidx = nothing
    end

    solAlg = AutoTsit5(Rodas5())
    sol = solve(prob, solAlg; saveat = tps, save_idxs = sidx, reltol = 1.0e-9, maxiters = 1e7)

    if sol.retcode != :Success
        println("Solving failed with the following parameters.")
        println(sol.retcode)
        println(ForwardDiff.value.(params))
        @assert sol.retcode == :Success
    end

    sol = hcat(sol.u...)'
    @assert ndims(sol) == 2
    @assert size(sol, 1) == length(tps)

    if pSTAT5
        @assert size(sol) == (length(tps), 3)
        # Summation of active species
        return vec(sol[:, 1] + 2 * (sol[:, 2] + sol[:, 3]))
    end

    return sol
end


" Calculate the Jacobian of the model and perform variance propagation. "
function runCkineVarProp(tps::Vector, params::Vector, sigma)::Vector
    checkInputs(tps, params)

    # Sigma is the covariance matrix of the input parameters
    function jacF(x)
        pp = vcat(params[1:24], x, params[30:end])
        return runCkine(tps, pp, pSTAT5 = true)
    end

    jac = zeros(5, length(tps))
    jacobian!(jac, jacF, params[25:29])

    # Just return the diagonal for the marginal variance
    return diag(transpose(jac) * sigma * jac)
end


include("fit.jl")
include("gprocess.jl")

export runCkine, runCkineVarProp, runCkineCost

end # module
