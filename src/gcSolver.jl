__precompile__()
module gcSolver

using OrdinaryDiffEq
using LinearAlgebra
using SteadyStateDiffEq
using DiffEqSensitivity
using Sundials
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

    sol = solve(prob, AutoTsit5(Rodas5()); saveat=tps, options...).u

    if length(tps) > 1
        sol = vcat(transpose.(sol)...)
    else
        sol = reshape(sol[1], (1, Nspecies))
    end

    return sol
end


function runCkineAS(tps::Vector{Float64}, params::Vector, reduce::Vector, data::Vector)
    @assert all(tps .>= 0.0)
    @assert length(tps) == length(data)
    u0 = solveAutocrine(params)

    prob = ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, AutoTsit5(Rodas5()); options...)

    dg = (out, u, p, t, i) -> out .= abs(data[i] - dot(u, reduce))
    du0, dp = adjoint_sensitivities_u0(sol, CVODE_BDF(), dg, tps)
    u0g = ForwardDiff.jacobian(solveAutocrine, params)

    return dp .+ (u0g' * du0)
end


export runCkine, runCkineAS

precompile(runCkine, (Array{Float64, 1}, Array{Float64, 1}))

end # module
