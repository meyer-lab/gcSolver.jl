const Nspecies = 47 # number of complexes in surface + endosome + free ligand
const halfL = 19 # number of complexes on surface alone
const internalFrac = 0.5 # Same as that used in TAM model
const recIDX = [1, 2, 3, 10, 17]
const recIDXint = [ii + halfL for ii in recIDX]
const ligIDX = [39, 40, 41]
const activeSpec = [8, 9, 15, 16, 19]
const pSTATidx = [43, 44, 45]

const Nparams = 36 # number of unknowns for the full model
const kfbnd = 0.60 # Assuming on rate of 10^7 M-1 sec-1
const internalV = 623.0 # Same as that used in TAM model

# p[1:4] is kfwd, k1rev, k2rev, k4rev
# p[5:8] is k5rev, k10rev, k11rev, k13rev
# p[9:12] is k14rev, k16rev, k17rev, k22rev
# p[13:14] is k23rev, k24rev

function dYdT(du, u, p, ILs)
    # IL2Ra, IL2Rb, gc, IL2_IL2Ra, IL2_IL2Rb, IL2_IL2Ra_IL2Rb, IL2_IL2Ra_gc, IL2_IL2Rb_gc, IL2_IL2Ra_IL2Rb_gc
    # IL15Ra, IL15_IL15Ra, IL15_IL2Rb, IL15_IL15Ra_IL2Rb, IL15_IL15Ra_gc, IL15_IL2Rb_gc, IL15_IL15Ra_IL2Rb_gc

    # IL2/15
    for i = 0:1
        pp = view(p, (1 + i * 6):(8 + i * 6))
        k12rev = pp[2] * pp[7] / pp[3] # Detailed balance
        k8rev = pp[6] * k12rev / pp[5] # Detailed balance
        k9rev = pp[6] * pp[7] / pp[4] # Detailed balance

        du[2] +=
            -kfbnd * u[2] * ILs[i + 1] + pp[3] * u[5 + 7 * i] - p[1] * u[2] * u[7 + 7 * i] + k9rev * u[9 + 7 * i] - p[1] * u[2] * u[4 + 7 * i] +
            pp[7] * u[6 + 7 * i]
        du[3] +=
            -p[1] * u[5 + 7 * i] * u[3] + p[5 + 7 * i] * u[8 + 7 * i] - p[1] * u[4 + 7 * i] * u[3] + p[4 + 7 * i] * u[7 + 7 * i] -
            p[1] * u[6 + 7 * i] * u[3] + p[6 + 7 * i] * u[9 + 7 * i]
        du[1 + 9 * i] =
            -kfbnd * u[1 + 9 * i] * ILs[i + 1] + pp[2] * u[4 + 7 * i] - p[1] * u[1 + 9 * i] * u[8 + 7 * i] + k8rev * u[9 + 7 * i] -
            p[1] * u[1 + 9 * i] * u[5 + 7 * i] + k12rev * u[6 + 7 * i]
        du[4 + 7 * i] =
            -p[1] * u[4 + 7 * i] * u[2] + pp[7] * u[6 + 7 * i] - p[1] * u[4 + 7 * i] * u[3] +
            pp[4] * u[7 + 7 * i] +
            kfbnd * ILs[i + 1] * u[1 + 9 * i] - pp[2] * u[4 + 7 * i]
        du[5 + 7 * i] =
            -p[1] * u[5 + 7 * i] * u[1 + 9 * i] + k12rev * u[6 + 7 * i] - p[1] * u[5 + 7 * i] * u[3] +
            p[5 + 7 * i] * u[8 + 7 * i] +
            kfbnd * ILs[i + 1] * u[2] - pp[3] * u[5 + 7 * i]
        du[6 + 7 * i] =
            -p[1] * u[6 + 7 * i] * u[3] + pp[6] * u[9 + 7 * i] + p[1] * u[4 + 7 * i] * u[2] - pp[7] * u[6 + 7 * i] +
            p[1] * u[5 + 7 * i] * u[1 + 9 * i] - k12rev * u[6 + 7 * i]
        du[9 + 7 * i] =
            p[1] * u[8 + 7 * i] * u[1 + 9 * i] - k8rev * u[9 + 7 * i] + p[1] * u[7 + 7 * i] * u[2] - k9rev * u[9 + 7 * i] +
            p[1] * u[6 + 7 * i] * u[3] - pp[6] * u[9 + 7 * i]
        du[7 + 7 * i] = -p[1] * u[7 + 7 * i] * u[2] + k9rev * u[9 + 7 * i] + p[1] * u[4 + 7 * i] * u[3] - pp[4] * u[7 + 7 * i]
        du[8 + 7 * i] = -p[1] * u[8 + 7 * i] * u[1 + 9 * i] + k8rev * u[9 + 7 * i] + p[1] * u[3] * u[5 + 7 * i] - pp[5] * u[8 + 7 * i]
    end

    du[3] += -p[1] * u[3] * u[18] + p[15] * u[19]
    du[17] = -kfbnd * u[17] * ILs[3] + p[14] * u[18]
    du[18] = kfbnd * u[17] * ILs[3] - p[14] * u[18] - p[1] * u[3] * u[18] + p[15] * u[19]
    du[19] = p[1] * u[3] * u[18] - p[15] * u[19]

    return nothing
end


function fullDeriv(du, u, p, t)
    fill!(du, 0.0)

    ILs = view(p, 1:3)
    pSurf = view(p, 4:18)
    pEndo = p[19]
    trafP = view(p, 20:30)

    # Calculate cell surface and endosomal reactions
    dYdT(du, u, pSurf, ILs)

    # Add STAT5 reactions
    activR = sum(u[activeSpec]) + internalFrac * sum(u[activeSpec .+ halfL])
    du[42] = p[36] * u[47] - p[31] * u[42] * activR # STAT5
    du[43] = p[31] * u[42] * activR - p[32] * u[43] # pSTAT5
    du[44] = 0.5 * p[32] * u[43] - p[33] * u[44] # pSTAT5d
    du[45] = p[33] * u[44] - p[34] * u[45] # pSTAT5nd
    du[46] = p[34] * u[45] - 0.5 * p[35] * u[46] # STAT5nd
    du[47] = p[35] * u[46] - p[36] * u[47] # STAT5n

    dYdT(view(du, (halfL + 1):(2 * halfL)), view(u, (halfL + 1):(2 * halfL)), pEndo * pSurf, view(u, ligIDX))

    # Handle endosomal ligand balance.
    # Must come before trafficking as we only calculate this based on reactions balance
    du[39] = -sum(view(du, (halfL + 4):(halfL + 9))) / internalV
    du[40] = -sum(view(du, (halfL + 11):(halfL + 16))) / internalV
    du[41] = -sum(view(du, (halfL + 18):(halfL + 19))) / internalV

    # Actually calculate the trafficking
    for ii = 1:halfL
        if ii ∈ activeSpec
            du[ii] -= u[ii] * (trafP[1] + trafP[2]) # Endo
            du[ii + halfL] += u[ii] * (trafP[1] + trafP[2]) / internalFrac - trafP[5] * u[ii + halfL] # Endo, deg
        else
            du[ii] += -u[ii] * trafP[1] + trafP[4] * (1.0 - trafP[3]) * u[ii + halfL] * internalFrac # Endo, recycling
            du[ii + halfL] += u[ii] * trafP[1] / internalFrac - trafP[4] * (1.0 - trafP[3]) * u[ii + halfL] - (trafP[5] * trafP[3]) * u[ii + halfL] # Endo, rec, deg
        end
    end

    # Expression: IL2Ra, IL2Rb, gc, IL15Ra, IL7Ra
    du[recIDX] .+= view(trafP, 6:10)

    # Degradation does lead to some clearance of ligand in the endosome
    du[ligIDX] .-= view(u, ligIDX) .* trafP[5]

    return nothing
end


" Takes in parameters and solves for steady state expression and initial species states "
function solveAutocrine(rIn::Vector)
    @assert all(rIn .>= 0.0)
    r = view(rIn, 20:30)
    @assert r[3] < 1.0

    # r is endo, activeEndo, sortF, kRec, kDeg, Rexpr*8
    y0 = zeros(eltype(r), Nspecies)

    # Add base amount of STAT5
    y0[42] = r[11]

    # Expand out trafficking terms
    kRec = r[4] * (1.0 - r[3])
    kDeg = r[5] * r[3]

    # Assuming no autocrine ligand, so can solve steady state
    # Add the species
    y0[recIDXint] = r[6:10] / kDeg / internalFrac
    y0[recIDX] = (r[6:10] + kRec * y0[recIDXint] * internalFrac) / r[1]

    return y0
end
