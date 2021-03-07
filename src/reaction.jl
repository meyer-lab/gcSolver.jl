using Catalyst

const Nspecies = 47 # number of complexes in surface + endosome + free ligand
const halfL = 11 # number of complexes on surface alone
const internalFrac = 0.5 # Same as that used in TAM model
const recIDX = [1, 2, 3]
const recIDXint = [ii + halfL for ii in recIDX]
const ligIDX = [39, 40]
const activeSpec = [7, 8, 11]
const pSTATidx = [43, 44, 45]

const Nparams = 36 # number of unknowns for the full model
const internalV = 623.0 # Same as that used in TAM model


# 0.6 is kfbnd assuming on rate of 10^7 M-1 sec-1
rxn = @reaction_network begin
    (0.6 * IL2, k1ᵣ), Ra ↔ IL2_Ra
    (0.6 * IL2, k2ᵣ), Rb ↔ IL2_Rb
    (kf, k4ᵣ), IL2_Ra + gc ↔ IL2_Ra_gc
    (kf, k5ᵣ), IL2_Rb + gc ↔ IL2_Rb_gc
    (kf, k1ᵣ * k10ᵣ * k11ᵣ / k2ᵣ / k5ᵣ), IL2_Rb_gc + Ra ↔ IL2_Ra_Rb_gc # k8r, detailed balance
    (kf, k10ᵣ * k11ᵣ / k4ᵣ), IL2_Ra_gc + Rb ↔ IL2_Ra_Rb_gc # k9r, detailed balance
    (kf, k10ᵣ), IL2_Ra_Rb + gc ↔ IL2_Ra_Rb_gc
    (kf, k11ᵣ), IL2_Ra + Rb ↔ IL2_Ra_Rb
    (kf, k1ᵣ * k11ᵣ / k2ᵣ), IL2_Rb + Ra ↔ IL2_Ra_Rb # k12r, detailed balance
    (0.6 * IL15, k14ᵣ), Rb ↔ IL15_Rb
    (kf, k17ᵣ), IL15_Rb + gc ↔ IL15_Rb_gc
end IL2 IL15 kf k1ᵣ k2ᵣ k4ᵣ k5ᵣ k10ᵣ k11ᵣ k14ᵣ k17ᵣ

println(rxn)

const rxnF = convert(ODESystem, rxn)


function fullDeriv(du, u, p, t)
    fill!(du, 0.0)

    pEndo = p[12]
    trafP = view(p, 20:30)

    # Calculate cell surface and endosomal reactions
    rxnF.f(du, u, p[1:11])
    endoIDX = (halfL + 1):(2 * halfL)
    rxnF.f(du[endoIDX], u[endoIDX], hcat(u[ligIDX], p[3:11] * pEndo))

    # Add STAT5 reactions
    activR = sum(u[activeSpec]) + internalFrac * sum(u[activeSpec .+ halfL])
    du[42] = p[36] * u[47] - p[31] * u[42] * activR # STAT5
    du[43] = p[31] * u[42] * activR - p[32] * u[43]^2 # pSTAT5
    du[44] = 0.5 * p[32] * u[43]^2 - p[33] * u[44] # pSTAT5d
    du[45] = p[33] * u[44] - p[34] * u[45] # pSTAT5nd
    du[46] = p[34] * u[45] - 0.5 * p[35] * u[46] # STAT5nd
    du[47] = p[35] * u[46] - p[36] * u[47] # STAT5n

    # Handle endosomal ligand balance.
    # Must come before trafficking as we only calculate this based on reactions balance
    du[39] = -sum(view(du, (halfL + 4):(halfL + 9))) / internalV
    du[40] = -sum(view(du, (halfL + 10):(halfL + 11))) / internalV

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
    du[ligIDX] .-= view(u, ligIDX) * trafP[5]

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
