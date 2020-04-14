using StaticArrays

const Nspecies = 47 # number of complexes in surface + endosome + free ligand
const halfL = 19 # number of complexes on surface alone
const internalFrac = 0.5 # Same as that used in TAM model
const recIDX = SVector(1, 2, 3, 10, 17)
const recIDXint = @SVector [ii + halfL for ii in recIDX]
const ligIDX = SVector(39, 40, 41)
const activeSpec = SVector(8, 9, 15, 16, 19)

const Nparams = 36 # number of unknowns for the full model
const Nlig = 3 # Number of ligands
const kfbnd = 0.60 # Assuming on rate of 10^7 M-1 sec-1
const internalV = 623.0 # Same as that used in TAM model

# p[1:4] is kfwd, k1rev, k2rev, k4rev
# p[5:8] is k5rev, k10rev, k11rev, k13rev
# p[9:12] is k14rev, k16rev, k17rev, k22rev
# p[13:14] is k23rev, k24rev
"""
    kfbnd = 0.60
    rd["IL2"], rd["IL15"], rd["IL7"], rd["IL9"], rd["IL4"], rd["IL21"], rd["kfwd"] = tuple(rxntfr[0:7])
    rd["surf.k1rev"] = kfbnd * 10.0  # 7
    rd["surf.k2rev"] = kfbnd * 144.0
    rd["surf.k4rev"], rd["surf.k5rev"] = rxntfr[7], rxntfr[8]  # 9 #10
    rd["surf.k10rev"] = 12.0 * rd["surf.k5rev"] / 1.5
    rd["surf.k11rev"] = 63.0 * rd["surf.k5rev"] / 1.5
    rd["surf.k13rev"] = kfbnd * 0.065
    rd["surf.k14rev"] = kfbnd * 438.0
    rd["surf.k16rev"] = rxntfr[9]
    rd["surf.k17rev"] = rxntfr[10]  # 16
    rd["surf.k22rev"] = rxntfr[11]
    rd["surf.k23rev"] = rxntfr[12]
    rd["surf.k25rev"] = kfbnd * 59.0
    rd["surf.k27rev"] = rxntfr[13]

"""
function getUnkVec():
    """Creates full vector of unknown values to be fit"""
    #kfwd, k4, k5, k16, k17, k22, k23, k27, endo, aendo, sort, krec, kdeg, k34, k35, k36, k37, k38, k39
    unkVecF = zeros(Float64, 1, 19)
    
    unkVecF[1] = .00125
    unkVecF[2:7] = 0.679 # Mean of previously used prior distribution (lognormal, mu=0, sd=0.5)
    unkVecF[8] = 1.0 # Previous prior
    unkVecF[9] = 0.1
    unkvecF[10] = 0.678
    unkvecF[11] = 0.1
    unk0.01
    asdfsdf0.13
    unasdfF[14:19] = 0.679, # Same as previous rationale
end
    

function fitParams(ILs, unkVec, recAbundances):
    """Takes in full unkvec and constructs it into full fit parameters vector"""
    kfbnd = 0.60
    paramvec = zeros(Float64, 1, Nparams)
    paramvec[1] = 
    
    paramvec
    
end


function yHatVec():
    """Constructs full vector of pSTAT means and variances to fit to"""
    #import data into Julia Vector - should be X by 2
end


function resids(x):
    #TODO add weights etc.
    fit = fitParams(x)
    ytrue = yHatVec()
    return runmodel(fitmat) - ytrue
end


function runFit():
    unkVecInit = getUnkVec
    optimize(resids, unkVecInit)
    return fit
end
    