module gcSolver

using OrdinaryDiffEq
using ForwardDiff
using StaticArrays

include("reaction.jl")

function fullDeriv(du, u, p, t)
    ILs, surface, endosome, trafP = fullParam(p)

    fullModel(du, u, surface, endosome, trafP, ILs)
end


function fullParam(rxntfR)
    surface = ones(eltype(rxntfR), 21)

    if length(rxntfR) == Nparams
        ILs = rxntfR[1:6]
        surface[[1, 4, 5]] = rxntfR[7:9]
        surface[2] = kfbnd * 10 # doi:10.1016/j.jmb.2004.04.038, 10 nM
        surface[3] = kfbnd * 144 # doi:10.1016/j.jmb.2004.04.038, 144 nM
        surface[7] = 63.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
        surface[8] = kfbnd * 0.065 # based on the multiple papers suggesting 30-100 pM
        surface[9] = kfbnd * 438 # doi:10.1038/ni.2449, 438 nM
        surface[10:13] = rxntfR[10:13]
        surface[14] = kfbnd * 59 # DOI:10.1111/j.1600-065X.2012.01160.x, 59 nM
        surface[[15, 17, 19, 21]] = view(rxntfR, 14:17)
        surface[16] = kfbnd * 0.1 # DOI:10.1073/pnas.89.12.5690, ~100 pM
        surface[18] = kfbnd * 1.0 # DOI: 10.1126/scisignal.aal1253 (human)
        surface[20] = kfbnd * 0.07 # DOI: 10.1126/scisignal.aal1253 (human)

        trafP = rxntfR[18:Nparams]
    else
        @assert length(rxntfR) == NIL2params
        ILs = zeros(eltype(rxntfR), Nlig)
        ILs[1] = rxntfR[1]
        surface[1:5] .= rxntfR[2:6]
        surface[6] = rxntfR[7]

        trafP = zeros(eltype(rxntfR), 13)
        trafP[1:5] = [0.08, 1.46, 0.18, 0.15, 0.017]
        trafP[6:8] = rxntfR[8:10]
    end

    surface[6] = 12.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038

    endosome = copy(surface)
    # all reverse rates are 5-fold higher in endosome
    endosome[2:21] *= 5.0

    @assert trafP[3] < 1.0
    @assert all(ILs .>= 0.0)
    @assert all(surface .>= 0.0)
    @assert all(endosome .>= 0.0)
    @assert all(trafP .>= 0.0)
    @assert length(ILs) == Nlig

    return ILs, surface, endosome, trafP
end


function runCkine(tps, params; alg=Rosenbrock23())
    _, _, _, trafP = fullParam(params)

    u0 = solveAutocrine(trafP)

    prob = ODEProblem(fullDeriv, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, alg; reltol=1.0e-6, abstol=1.0e-6, isoutofdomain=(u, p, t) -> any(x -> x < 0.0, u))

    return sol(tps).u
end


function runCkineS(tps, params)
    return ForwardDiff.gradient(pp -> runCkine(tps, pp), params)::Vector{Float64}
end


export runCkine
export runCkineS


end # module
