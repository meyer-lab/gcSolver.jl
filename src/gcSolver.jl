
__precompile__()
module gcSolver

using OrdinaryDiffEq
using LinearAlgebra
using SteadyStateDiffEq
using DiffEqSensitivity

include("reaction.jl")

function fullDeriv(du, u, p, t)
    fill!(du, 0.0)

    fullModel(du, u, view(p, 7:27), view(p, 28:48), trafP(p), view(p, 1:6))
end


const solTol = 1.0e-9


function domainDef(u, p, t)
    return any(x -> x < -solTol, u)
end


function runCkine(tps::Vector{Float64}, params::Vector)::Matrix
    @assert all(params .>= 0.0)
    @assert all(tps .>= 0.0)

    u0 = solveAutocrine(params)

    prob = ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, AutoTsit5(Rodas5()); reltol = solTol, abstol = solTol, isoutofdomain = domainDef)
    solut = sol(tps).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, Nspecies))
    end

    return solut
end



function runCkineFS(tp::Real, params::Vector)
    @assert all(params .>= 0.0)
    @assert tp > 0.0

    u0 = solveAutocrine(params)

    prob = ODEForwardSensitivityProblem(fullDeriv, u0, (0.0, tp), params)

    sol = solve(prob, Rosenbrock32(autodiff=false); reltol = solTol, abstol = solTol, isoutofdomain = domainDef)

    return extract_local_sensitivities(sol, tp)
end


# TODO: Not finished.
function runCkineAS(tps::Vector{Float64}, params::Vector, reduce::Vector, data::Vector)
    @assert all(params .>= 0.0)
    @assert tp > 0.0

    dg = (out, u, i) -> out .= data[i] .- dot(u, reduce)

    u0 = solveAutocrine(params)

    prob = ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, AutoTsit5(Rodas5()); reltol = solTol, abstol = solTol, isoutofdomain = domainDef)

    res = adjoint_sensitivities(sol, AutoTsit5(Rodas5()), dg, tps, abstol=solTol, reltol=solTol, iabstol=solTol, ireltol=solTol)

    return res
end



function runCkineSS(params::Vector)
    @assert all(params .>= 0.0)

    u0 = solveAutocrine(params)

    probInit = SteadyStateProblem(fullDeriv, u0, params)

    solInit = solve(probInit, DynamicSS(AutoTsit5(Rodas5())); reltol = solTol, abstol = solTol, isoutofdomain = domainDef)

    return solInit
end


export runCkine, runCkineSS, runCkineFS

precompile(runCkine, (Array{Float64, 1}, Array{Float64, 1}))
precompile(runCkineSS, (Array{Float64, 1},))

end # module
