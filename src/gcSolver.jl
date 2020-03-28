__precompile__()
module gcSolver

using OrdinaryDiffEq
using LinearAlgebra
using SteadyStateDiffEq
using ForwardDiff

include("reaction.jl")

function fullDeriv(du, u, p, t)
    fill!(du, 0.0)

    fullModel(du, u, view(p, 7:27), view(p, 28:48), trafP(p), view(p, 1:6))
end


const solTol = 1.0e-9

function domainDef(u, p, t)
    return any(x -> x < -solTol, u)
end

const options = Dict([:reltol => solTol, :abstol => solTol, :isoutofdomain => domainDef])


function runCkine(tps::Vector{Float64}, params::Vector)::Matrix
    @assert all(tps .>= 0.0)
    u0 = solveAutocrine(params)

    prob = ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), params)

    alg = AutoTsit5(Rodas5(autodiff = eltype(params) == Float64))
    sol = solve(prob, alg; saveat = tps, options...).u

    if length(tps) > 1
        sol = vcat(transpose.(sol)...)
    else
        sol = reshape(sol[1], (1, Nspecies))
    end

    return sol
end


" Converts the ODE solution to a predicted amount of pSTAT. "
function runCkinePSTAT(tps::Vector, params::Vector)
    # TODO: Add in the sigmoidal relationship.
    retval = runCkine(tps::Vector{Float64}, params::Vector)

    # Summation of active species
    pSTAT = sum(retval[:, SVector(8, 9, 15, 16, 19, 22, 25, 28)], dims=2)

    @assert length(pSTAT) == length(tps)
    return vec(pSTAT)
end


" Calculate the Jacobian of the model and perform variance propagation. "
function runCkineVarPorp(tps::Vector, params::Vector, sigma)::Matrix
    function jacF(x)
        return runCkinePSTAT(tps, x)
    end

    jac = zeros(length(params), length(tps))
    ForwardDiff.jacobian!(jac, jacF, params)

    return transpose(jac) * sigma * jac
end


include("normFlows.jl")

export runCkine

precompile(runCkine, (Array{Float64, 1}, Array{Float64, 1}))

end # module
