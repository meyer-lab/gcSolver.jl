module gcSolver

using OrdinaryDiffEq
import LinearAlgebra: diag, norm
import ForwardDiff
using Optim
using Statistics
import ModelingToolkit

include("reaction.jl")

const solTol = 1.0e-5

function domainDef(u, p, t)
    return any(x -> x < 0.0, u)
end


" Check that we haven't been provided anything ridiculous by the user. "
function checkInputs(tps::Vector{Float64}, params::Vector)
    @assert all(tps .>= 0.0)
    @assert params[22] < 1.0
    @assert length(params) == Nparams
    @assert all(params .>= 0.0)
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
        sidx = @SVector [43, 44, 45]
    else
        sidx = nothing
    end

    sol = solve(prob, AutoTsit5(KenCarp5(), stiffalgfirst = true); saveat = tps, reltol = solTol, save_idxs = sidx, isoutofdomain = domainDef).u

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
        return vec(sol[:, 1] + 2 * (sol[:, 2] + sol[:, 3]))
    end

    return sol
end


" Calculate the Jacobian of the model and perform variance propagation. "
function runCkineVarProp(tps::Vector, params::Vector, sigma)::Vector
    checkInputs(tps, params)

    # Sigma is the covariance matrix of the input parameters
    function jacF(x)
        pp = vcat(params[1:27], x, params[33:end])
        return runCkine(tps, pp, pSTAT5 = true)
    end

    jac = zeros(5, length(tps))
    ForwardDiff.jacobian!(jac, jacF, params[28:32])

    # Just return the diagonal for the marginal variance
    return diag(transpose(jac) * sigma * jac)
end


" Calculate the Hessian of the model. "
function runCkineHessian(tps::Vector, params::Vector)::Array
    checkInputs(tps, params)

    # Sigma is the covariance matrix of the input parameters
    function hF(tt::Float64, x)::Real
        pp = vcat(params[1:27], x, params[33:end])
        return runCkine([tt], pp, pSTAT5 = true)[1]
    end

    H = zeros(5, 5, length(tps))

    for ii = 1:length(tps)
        ForwardDiff.hessian!(view(H, :, :, ii), (x) -> hF(tps[ii], x), params[28:32])
    end

    return H
end


include("fit.jl")

export runCkine, runCkineVarProp

end # module
