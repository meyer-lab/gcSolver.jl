using CSV
using Memoize

const dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")


"""Creates full vector of unknown values to be fit"""
function getUnkVec()
    #kfwd, k4, k5, k16, k17, k22, k23, k27, endo, aendo, sort, krec, kdeg, k34, k35, k36, k37, k38, k39
    unkVecF = zeros(19)

    unkVecF[1] = 0.00125 # means of prior distributions from gc-cytokines paper
    unkVecF[2:7] .= 0.679
    unkVecF[8] = 1.0
    unkVecF[9] = 0.1
    unkVecF[10] = 0.678
    unkVecF[11] = 0.1
    unkVecF[12] = 0.01
    unkVecF[13] = 0.13
    unkVecF[14:19] .= 0.679 # pSTAT Rates

    return unkVecF
end


"""Takes in full unkvec and constructs it into full fit parameters vector - TODO move this"""
function fitParams(ILs, unkVec, recAbundances)
    kfbnd = 0.60
    paramvec = zeros(Float64, 1, Nparams)
    paramvec[1:3] = ILs
    paramvec[4] = unkVec[1] #kfwd
    paramvec[5] = kfbnd * 10.0 #k1rev
    paramvec[6] = kfbnd * 144.0 #k2rev
    paramvec[7] = unkVec[2] #k4rev
    paramvec[8] = unkVec[3] #k5rev
    paramvec[9] = 12.0 * paramvec[8] / 1.5 #k10
    paramvec[10] = 63.0 * paramvec[8] / 1.5 #k11
    paramvec[11] = kfbnd * 0.065 #k13
    paramvec[12] = kfbnd * 438.0 #k14
    paramvec[13] = unkVec[4] #k16
    paramvec[14] = unkVec[5] #k16
    paramvec[15] = unkVec[6] #k17
    paramvec[16] = unkVec[7] #k22
    paramvec[17] = kfbnd * 59.0 #k23
    paramvec[18] = unkVec[8] #k25
    paramvec[19] = 5.0 #endoadjust
    paramvec[20:24] = unkVec[9:13]
    paramvec[25:29] = receptorExp(recAbundances, unkvec[9], unkvec[11], unkvec[12], unkvec[13])
    paramvec[30] = 1.0 #TODO add initial STAT
    paramVec[31:36] = unkVec[14:19]

end


""" Uses receptor abundance (from flow) and trafficking rates to calculate receptor expression rate at steady state. """
function receptor_expression(receptor_abundance, endo, kRec, sortF, kDeg)
    rec_ex = (receptor_abundance * endo) / (1.0 + ((kRec * (1.0 - sortF)) / (kDeg * sortF)))
    return rec_ex
end


"""Constructs full vector of pSTAT means and variances to fit to, and returns expression levels for use with fitparams"""
@memoize function getyVec()
    #import data into Julia Vector - should be X by 2

    df = CSV.read(joinpath(dataDir, "VarianceData.csv"), copycols = true)
    sort!(df, (:Date, :Ligand, :Cell, :Dose, :Time))
    yVec = df.Mean #add in variance later
    cellVec = df.Cell
    tpsVec = df.Time * 60.0
    ligs = df.Ligand
    doses = df.Dose

    ligVec = zeros(Float64, size(df)[1], 3)
    expVec = zeros(Float64, size(df)[1], 5)
    expDict = getExpression()

    for i = 1:size(df)[1]
        if ligs[i] == "IL2"
            ligVec[i, 1] = doses[i]
        else
            ligVec[i, 2] = doses[i]
        end
        expVec[i, 1:5] = expDict(cellVec[i])
    end

    return yVec, tpsVec, expVec, ligVec
end


"""Gets expression vector for each cell type and puts it into dictionary"""
function getExpression()
    #CSV.read...
    #dict = {Treg, TregNaive,....}
    #entries = ...
    return recDict
end


"""Calculates squared error for a given unkVec"""
function resids(x::Vector{T}) where T
    #TODO add weights etc.
    ytrue, tps, expVec, ligVec = getyVec()
    yhat = similar(ytrue, T)
    for i = 1:size(tps)[1]
        vec = fitParams(ligVec[i, 1:3], x, expVec[i, 1:5])
        println(tps[i])
        yhat[i] = runCkinePSTAT(tps[i], vec)
    end

    return norm(yhat .- ytrue)
end


""" Gets inital unkowns, optimizes them, and returns parameters of best fit"""
function runFit(; itern = 1E6)
    unkVecInit = getUnkVec()
    fit = optimize(resids, unkVecInit, LBFGS(), Optim.Options(iterations = itern, show_trace = true))

    return fit.minimizer
end
