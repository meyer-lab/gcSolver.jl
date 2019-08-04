module gcSolver

using OrdinaryDiffEq
using DiffEqSensitivity

include("reaction.jl")

function fullDeriv(du, u, p, t)
    ILs, surface, endosome, trafP = fullParam(p)

    fullModel(du, u, surface, endosome, trafP, ILs)
end


function IL2Deriv(du::Vector, u::Vector, p::Vector, t)
    ILs, surface, endosome, trafP = IL2param(p)

    fullModel(du, u, surface, endosome, trafP, ILs)
end


function fullParam(rxntfR)
    ILs = rxntfR[1:6]
    surface = Vector{eltype(rxntfR)}(undef, 21)
    surface[[1, 4, 5]] = view(rxntfR, 7:9)
    surface[2] = kfbnd * 10 # doi:10.1016/j.jmb.2004.04.038, 10 nM
    surface[3] = kfbnd * 144 # doi:10.1016/j.jmb.2004.04.038, 144 nM
    surface[6] = 12.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
    surface[7] = 63.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038
    surface[8] = kfbnd * 0.065 # based on the multiple papers suggesting 30-100 pM
    surface[9] = kfbnd * 438 # doi:10.1038/ni.2449, 438 nM
    surface[10:13] = rxntfR[10:13]
    surface[14] = kfbnd * 59 # DOI:10.1111/j.1600-065X.2012.01160.x, 59 nM
    surface[[15, 17, 19, 21]] = view(rxntfR, 14:17)
    surface[16] = kfbnd * 0.1 # DOI:10.1073/pnas.89.12.5690, ~100 pM
    surface[18] = kfbnd * 1.0 # DOI: 10.1126/scisignal.aal1253 (human)
    surface[20] = kfbnd * 0.07 # DOI: 10.1126/scisignal.aal1253 (human)

    trafP = view(rxntfR, 18:Nparams)

    endosome = copy(surface)
    endosome[2:21] *= 5.0

    @assert trafP[3] < 1.0

    # all reverse rates are same in the endosome
    return ILs, surface, endosome, trafP
end


function IL2param(rxntfR)
    ILs = zeros(eltype(rxntfR), Nlig)
    ILs[1] = rxntfR[1]
    surface = ones(eltype(rxntfR), 21)
    surface[1:6] .= rxntfR[2:7]

    # These are probably measured in the literature
    surface[6] = 12.0 * surface[5] / 1.5 # doi:10.1016/j.jmb.2004.04.038

    trafP = zeros(eltype(rxntfR), 13)
    trafP[1:5] = [0.08, 1.46, 0.18, 0.15, 0.017]
    trafP[6:8] = rxntfR[8:10]

    endosome = copy(surface)
    endosome[2:21] *= 5.0

    @assert trafP[3] < 1.0

    return ILs, surface, endosome, trafP
end


function runCkine(tps, params, IL2case)
    if IL2case
        _, _, _, trafP = IL2param(params)
        f = IL2Deriv
    else
        _, _, _, trafP = fullParam(params)
        f = fullDeriv
    end

    u0 = solveAutocrine(trafP)

    prob = ODEProblem(f, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, Rosenbrock23())

    return sol(tps).u
end


export runCkine


end # module
