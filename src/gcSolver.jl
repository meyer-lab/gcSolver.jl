
__precompile__()
module gcSolver

using OrdinaryDiffEq
using ForwardDiff
using LinearAlgebra
using SteadyStateDiffEq

include("reaction.jl")

function fullDeriv(du, u, (p, surface, endosome, trafP, ILs), t)
    fill!(du, 0.0)
    fullParam!(p, surface, endosome, trafP, ILs)

    fullModel(du, u, surface, endosome, trafP, ILs)
end


function fullParam!(rxntfR::Vector, surface, endosome, trafP, ILs)
    fill!(surface, 1.0)
    fill!(endosome, 1.0)
    fill!(trafP, 0.0)
    fill!(ILs, 0.0)

    ILs[:] .= rxntfR[1:6]
    surface[SVector(1, 4, 5, 10, 11, 12, 13)] = view(rxntfR, 7:13)
    surface[2] = kfbnd * 10 # doi:10.1016/j.jmb.2004.04.038, 10 nM
    surface[3] = kfbnd * 144 # doi:10.1016/j.jmb.2004.04.038, 144 nM
    surface[6] = 12.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
    surface[7] = 63.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
    surface[8] = kfbnd * 0.065 # based on the multiple papers suggesting 30-100 pM
    surface[9] = kfbnd * 438 # doi:10.1038/ni.2449, 438 nM
    surface[14] = kfbnd * 59 # DOI:10.1111/j.1600-065X.2012.01160.x, 59 nM
    surface[SVector(15, 17, 19, 21)] = view(rxntfR, 14:17)
    surface[16] = kfbnd * 0.1 # DOI:10.1073/pnas.89.12.5690, ~100 pM
    surface[18] = kfbnd * 1.0 # DOI: 10.1126/scisignal.aal1253 (human)
    surface[20] = kfbnd * 0.07 # DOI: 10.1126/scisignal.aal1253 (human)

    # all reverse rates are 5-fold higher in endosome
    endosome[:] .= surface
    endosome[2:21] *= 5.0

    trafP[:] = rxntfR[18:Nparams]

    @assert trafP[3] < 1.0

    return ILs, surface, endosome, trafP
end


function domainDef(u, p, t)
    return any(x -> x < 0.0, u)
end


function runCkine(tps::Vector{Float64}, params::Vector)::Matrix{Float64}
    @assert all(params .>= 0.0)
    @assert all(tps .>= 0.0)

    # Allocate temporaries
    surface = @MVector ones(eltype(params), 21)
    endosome = @MVector ones(eltype(params), 21)
    ILs = @MVector zeros(eltype(params), Nlig)
    trafP = @MVector zeros(eltype(params), 13)

    fullParam!(params, surface, endosome, trafP, ILs)
    u0 = solveAutocrine(trafP)

    prob = ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), (params, surface, endosome, trafP, ILs))

    sol = solve(prob, Rodas5(); reltol = 1.0e-8, abstol = 1.0e-8, isoutofdomain = domainDef)
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

    # Allocate temporaries
    surface = @MVector ones(eltype(params), 21)
    endosome = @MVector ones(eltype(params), 21)
    ILs = @MVector zeros(eltype(params), Nlig)
    trafP = @MVector zeros(eltype(params), 13)

    fullParam!(params, surface, endosome, trafP, ILs)
    u0 = solveAutocrine(trafP)

    probInit = SteadyStateProblem(fullDeriv, u0, (params, surface, endosome, trafP, ILs))

    solInit = solve(probInit, DynamicSS(Rodas5()); isoutofdomain = domainDef)

    return solInit
end

export runCkine, runCkineSS

precompile(runCkine, (Array{Float64, 1}, Array{Float64, 1}))
precompile(runCkineSS, (Array{Float64, 1},))

end # module
