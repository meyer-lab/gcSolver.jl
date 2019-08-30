module gcSolver

using OrdinaryDiffEq
using ForwardDiff

include("reaction.jl")

function fullDeriv(du::Vector, u::Vector, p::Vector, t)
    fill!(du, 0.0)
    ILs, surface, endosome, trafP = fullParam(p)

    fullModel(du, u, surface, endosome, trafP, ILs)
end


function fullParam(rxntfR::Vector)
    surface = ones(eltype(rxntfR), 21)
    ILs = zeros(eltype(rxntfR), Nlig)
    trafP = zeros(eltype(rxntfR), 13)

    if length(rxntfR) == Nparams
        ILs[:] .= rxntfR[1:6]
        surface[[1, 4, 5, 10, 11, 12, 13]] = rxntfR[7:13]
        surface[2] = kfbnd * 10 # doi:10.1016/j.jmb.2004.04.038, 10 nM
        surface[3] = kfbnd * 144 # doi:10.1016/j.jmb.2004.04.038, 144 nM
        surface[6] = 12.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
        surface[7] = 63.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
        surface[8] = kfbnd * 0.065 # based on the multiple papers suggesting 30-100 pM
        surface[9] = kfbnd * 438 # doi:10.1038/ni.2449, 438 nM
        surface[14] = kfbnd * 59 # DOI:10.1111/j.1600-065X.2012.01160.x, 59 nM
        surface[[15, 17, 19, 21]] = view(rxntfR, 14:17)
        surface[16] = kfbnd * 0.1 # DOI:10.1073/pnas.89.12.5690, ~100 pM
        surface[18] = kfbnd * 1.0 # DOI: 10.1126/scisignal.aal1253 (human)
        surface[20] = kfbnd * 0.07 # DOI: 10.1126/scisignal.aal1253 (human)

        endosome = copy(surface)
        # all reverse rates are 5-fold higher in endosome
        endosome[2:21] *= 5.0

        trafP[:] = rxntfR[18:Nparams]
    else
        @assert length(rxntfR) == NIL2params
        ILs[1] = rxntfR[1]
        surface[1:5] .= rxntfR[2:6]
        surface[6] = 12.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
        surface[7] = rxntfR[7]

        endosome = copy(surface)
        # all reverse rates are 5-fold higher in endosome
        endosome[2:21] *= 5.0
        endosome[2:5] .= rxntfR[11:14]
        endosome[6] = 12.0 * endosome[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
        endosome[7] = rxntfR[15]

        trafP[1:5] = [0.08, 1.46, 0.18, 0.15, 0.017]
        trafP[6:8] = rxntfR[8:10]
    end

    @assert trafP[3] < 1.0

    return ILs, surface, endosome, trafP
end


function runCkine(tps::Array{Float64,1}, params::Vector)
    @assert all(params .>= 0.0)
    @assert all(tps .>= 0.0)
    _, _, _, trafP = fullParam(params)

    u0 = solveAutocrine(trafP)

    prob = ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, Rosenbrock23(); reltol=1.0e-6, abstol=1.0e-6, isoutofdomain=(u, p, t) -> any(x -> x < 0.0, u))

    solut = sol(tps).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, Nspecies))
    end

    return solut
end


function runCkineS(tps::Array{Float64,1}, params::Array{Float64,1})
    return ForwardDiff.jacobian(pp -> runCkine(tps, pp), params)
end


export runCkine, runCkineS

end # module
