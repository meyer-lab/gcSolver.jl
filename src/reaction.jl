using StaticArrays

const Nspecies = 25 # number of complexes in surface + endosome + free ligand
const halfL = 9 # number of complexes on surface alone
const internalFrac = 0.5 # Same as that used in TAM model
const recIDX = SVector(1, 2, 3)
const recIDXint = @SVector [ii + halfL for ii in recIDX]
const ligIDX = 19
const activeSpec = SVector(8, 9)

const Nparams = 24 # number of unknowns for the full model
const Nlig = 1 # Number of ligands
const kfbnd = 0.60 # Assuming on rate of 10^7 M-1 sec-1
const internalV = 623.0 # Same as that used in TAM model

# p[1:8] is kfwd, k1rev, k2rev, k4rev, k5rev, k10rev, k11rev, k13rev

function dYdT(du, u, p, ILs)
    # IL2Ra, IL2Rb, gc, IL2_IL2Ra, IL2_IL2Rb, IL2_IL2Ra_IL2Rb, IL2_IL2Ra_gc, IL2_IL2Rb_gc, IL2_IL2Ra_IL2Rb_gc

    # IL2
    k12rev = p[2] * p[7] / p[3] # Detailed balance
    k8rev = p[6] * k12rev / p[5] # Detailed balance
    k9rev = p[6] * p[7] / p[4] # Detailed balance

    du[2] += -kfbnd * u[2] * ILs[1] + p[3] * u[5] - p[1] * u[2] * u[7] + k9rev * u[9] - p[1] * u[2] * u[4] + p[7] * u[6]
    du[3] += -p[1] * u[5] * u[3] + p[5] * u[8] - p[1] * u[4] * u[3] + p[4] * u[7] - p[1] * u[6] * u[3] + p[6] * u[9]
    du[1] = -kfbnd * u[1] * ILs[1] + p[2] * u[4] - p[1] * u[1] * u[8] + k8rev * u[9] - p[1] * u[1] * u[5] + k12rev * u[6]
    du[4] = -p[1] * u[4] * u[2] + p[7] * u[6] - p[1] * u[4] * u[3] + p[4] * u[7] + kfbnd * ILs[1] * u[1] - p[2] * u[4]
    du[5] = -p[1] * u[5] * u[1] + k12rev * u[6] - p[1] * u[5] * u[3] + p[5] * u[8] + kfbnd * ILs[1] * u[2] - p[3] * u[5]
    du[6] = -p[1] * u[6] * u[3] + p[6] * u[9] + p[1] * u[4] * u[2] - p[7] * u[6] + p[1] * u[5] * u[1] - k12rev * u[6]
    du[9] = p[1] * u[8] * u[1] - k8rev * u[9] + p[1] * u[7] * u[2] - k9rev * u[9] + p[1] * u[6] * u[3] - p[6] * u[9]
    du[7] = -p[1] * u[7] * u[2] + k9rev * u[9] + p[1] * u[4] * u[3] - p[4] * u[7]
    du[8] = -p[1] * u[8] * u[1] + k8rev * u[9] + p[1] * u[3] * u[5] - p[5] * u[8]

    return nothing
end


function fullDeriv(du, u, p, t)
    fill!(du, 0.0)

    ILs = view(p, 1)
    pSurf = view(p, 2:8)
    pEndo = p[9]
    trafP = view(p, 10:18)
    STATp = view(p, 19:24)

    # Calculate cell surface and endosomal reactions
    dYdT(du, u, pSurf, ILs)

    # Add STAT5 reactions
    activR = sum(u[activeSpec]) + internalFrac * sum(u[activeSpec .+ halfL])
    du[20] = STATp[6] * u[25] - STATp[1] * u[20] * activR # STAT5
    du[21] = STATp[1] * u[20] * activR - STATp[2] * u[21] # pSTAT5
    du[22] = 0.5 * STATp[2] * u[21] - STATp[3] * u[22] # pSTAT5d
    du[23] = STATp[3] * u[22] - STATp[4] * u[23] # pSTAT5nd
    du[24] = STATp[4] * u[23] - 0.5 * STATp[5] * u[24] # STAT5nd
    du[25] = STATp[5] * u[24] - STATp[6] * u[25] # STAT5n

    dYdT(view(du, (halfL + 1):(2 * halfL)), view(u, (halfL + 1):(2 * halfL)), pEndo * pSurf, view(u, ligIDX))

    # Handle endosomal ligand balance.
    # Must come before trafficking as we only calculate this based on reactions balance
    du[ligIDX] = -sum(view(du, (halfL + 4):(halfL + 9))) / internalV

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

    # Expression: IL2Ra, IL2Rb, gc
    du[recIDX] .+= view(trafP, 6:8)

    # Degradation does lead to some clearance of ligand in the endosome
    du[ligIDX] -= u[ligIDX] * trafP[5]

    return nothing
end


" Takes in parameters and solves for steady state expression and initial species states. "
function solveAutocrine(rIn::Vector)
    @assert all(rIn .>= 0.0)
    r = view(rIn, 10:18)
    @assert r[3] < 1.0

    # r is endo, activeEndo, sortF, kRec, kDeg, Rexpr*8
    y0 = zeros(eltype(r), Nspecies)

    # Add base amount of STAT5
    y0[20] = r[9]

    # Expand out trafficking terms
    kRec = r[4] * (1.0 - r[3])
    kDeg = r[5] * r[3]

    # Assuming no autocrine ligand, so can solve steady state
    # Add the species
    y0[recIDXint] = r[6:8] / kDeg / internalFrac
    y0[recIDX] = (r[6:8] + kRec * y0[recIDXint] * internalFrac) / r[1]

    return y0
end
