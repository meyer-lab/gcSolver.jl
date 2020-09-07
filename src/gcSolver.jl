module gcSolver

using OrdinaryDiffEq
import LinearAlgebra: diag, norm, dot
import ForwardDiff
import ForwardDiff: value, Dual, partials, jacobian!, jacobian
using Optim
using Statistics
import ModelingToolkit: modelingtoolkitize
using DiffEqSensitivity
using Gadfly
gdf = Gadfly
using Plots
plt = Plots
import CSV
using DataFrames
import StatsBase: indicatormat
using StatsFuns

include("reaction.jl")
include("dataImport.jl")

const solTol = 1.0e-9
const solAlg = AutoTsit5(KenCarp4(), stiffalgfirst = true)


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


function runCkineCost(tps::Vector{Float64}, params::Vector{Dual{T, V, N}}, dataa) where {T, V, N}
    prob = runCkineSetup(tps, value.(params))
    sAlg = InterpolatingAdjoint()

    sol = solve(prob, solAlg; reltol = 1.0e-9)
    dg = (out, u, p, t, i) -> out .= dataa[i] - (u[pSTATidx[1]] + 2 * (u[pSTATidx[2]] + u[pSTATidx[3]]))
    du0, dp = adjoint_sensitivities(sol, KenCarp4(autodiff = false), dg, tps; sensealg = sAlg)

    # Convert du0 into parameter effects
    J = jacobian(solveAutocrine, value.(params))
    dOverall = vec(du0' * J + dp)

    sol = sol(tps).u
    if length(tps) > 1
        sol = vcat(transpose.(sol)...)
    else
        sol = reshape(sol[1], (1, Nspecies))
    end

    soll = vec(sol[:, pSTATidx[1]] + 2 * (sol[:, pSTATidx[2]] + sol[:, pSTATidx[3]]))
    cost = norm(soll - dataa)

    part = partials(params[1]) * dOverall[1]
    for ii = 2:length(dOverall)
        part += partials(params[ii]) * dOverall[ii]
    end

    return Dual{T, V, N}(cost, part)
end


" Actually run the gc ODE model. "
function runCkine(tps::Vector{Float64}, params; pSTAT5 = false)
    if params isa Vector
        prob = runCkineSetup(tps, params)
    else
        prob = params
    end

    sol = solve(prob, solAlg; saveat = tps, reltol = 1.0e-9).u

    if length(tps) > 1
        sol = vcat(transpose.(sol)...)
    else
        sol = reshape(sol[1], (1, Nspecies))
    end

    if length(tps) > size(sol, 1)
        println("Solving failed with the following parameters.")
        println(ForwardDiff.value.(params))
        @assert length(tps) == size(sol, 1)
    end

    if pSTAT5
        # Summation of active species
        return vec(sol[:, pSTATidx[1]] + 2 * (sol[:, pSTATidx[2]] + sol[:, pSTATidx[3]]))
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
