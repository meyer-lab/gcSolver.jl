
__precompile__()
module gcSolver

using OrdinaryDiffEq
using ForwardDiff
using LinearAlgebra
using SteadyStateDiffEq

include("reaction.jl")

function fullDeriv(du, u, p, t)
    fill!(du, 0.0)

    fullModel(du, u, view(p, 7:27), view(p, 28:48), trafP(p), view(p, 1:6))
end


function domainDef(u, p, t)
    return any(x -> x < 0.0, u)
end


const solTol = 1.0e-8

function runCkine(tps::Vector{Float64}, params::Vector)::Matrix{Float64}
    @assert all(params .>= 0.0)
    @assert all(tps .>= 0.0)

    u0 = solveAutocrine(params)

    prob = ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, ABDF2(); reltol = solTol, abstol = solTol, isoutofdomain = domainDef)
    solut = sol(tps).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, Nspecies))
    end

    return solut
end


function runCkineSS(params::Vector)
    @assert all(params .>= 0.0)

    u0 = solveAutocrine(params)

    probInit = SteadyStateProblem(fullDeriv, u0, params)

    solInit = solve(probInit, DynamicSS(ABDF2()); reltol = solTol, abstol = solTol, isoutofdomain = domainDef)

    return solInit
end

export runCkine, runCkineSS

precompile(runCkine, (Array{Float64, 1}, Array{Float64, 1}))
precompile(runCkineSS, (Array{Float64, 1},))

end # module
