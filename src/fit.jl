using CSV
using Memoize

const dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")


"""Creates full vector of unknown values to be fit"""
function getUnkVec()
    #kfwd, k4, k5, k16, k17, k22, k23, k27, endo, aendo, sort, krec, kdeg, k34, k35, k36, k37, k38, k39
    unkVecF = zeros(20)

    unkVecF[1] = 0.00125 # means of prior distributions from gc-cytokines paper
    unkVecF[2:7] .= 0.679
    unkVecF[8] = 1.0
    unkVecF[9] = 0.1
    unkVecF[10] = 0.678
    unkVecF[11] = 0.1
    unkVecF[12] = 0.01
    unkVecF[13] = 0.13
    unkVecF[14] = 30000.0
    unkVecF[15:20] .= 0.679 # pSTAT Rates

    return unkVecF
end


"""Takes in full unkvec and constructs it into full fit parameters vector - TODO move this"""
function fitParams(ILs, unkVec::Vector{T}, recAbundances) where {T}
    kfbnd = 0.60
    paramvec = zeros(T, Nparams)
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
    paramvec[13:16] = unkVec[4:7] #k16, k16, k17, k22
    paramvec[17] = kfbnd * 59.0 #k23
    paramvec[18] = unkVec[8] #k25
    paramvec[19] = 5.0 #endoadjust
    paramvec[20:24] = unkVec[9:13]
    paramvec[25:29] = receptor_expression(recAbundances, unkVec[9], unkVec[11], unkVec[12], unkVec[13])
    paramvec[30] = unkVec[14]
    paramvec[31:36] = unkVec[15:20]

    return paramvec
end


""" Uses receptor abundance (from flow) and trafficking rates to calculate receptor expression rate at steady state. """
function receptor_expression(receptor_abundance, endo, kRec, sortF, kDeg)
    rec_ex = (receptor_abundance * endo) / (1.0 + ((kRec * (1.0 - sortF)) / (kDeg * sortF)))
    return rec_ex
end


"""Constructs full vector of pSTAT means and variances to fit to, and returns expression levels for use with fitparams"""
@memoize function getyVec()
    #import data into Julia Vector - should be X by 1 until variance is added

    df = CSV.read(joinpath(dataDir, "MomentFitData.csv"), copycols = true)
    sort!(df, (:Date, :Ligand, :Cell, :Dose, :Time))
    yVec = df.Mean #add in variance later
    cellVec = df.Cell
    tpsVec = vec(df.Time * 60.0)
    ligs = df.Ligand
    doses = df.Dose

    ligVec = zeros(size(df)[1], 3)
    expVec = zeros(size(df)[1], 5)
    exprDF = getExpression()

    for i = 1:size(df)[1]
        if ligs[i] == "IL2"
            ligVec[i, 1] = doses[i]
        else
            ligVec[i, 2] = doses[i]
        end
        expVec[i, 1:5] = exprDF[!, Symbol(cellVec[i])]
    end

    return yVec, tpsVec, expVec, ligVec
end


"""Gets expression vector for each cell type and puts it into dictionary"""
function getExpression()
    recDF = CSV.read(joinpath(dataDir, "FakeExpressionData.csv"), copycols = true)
    return recDF
end


"""Calculates squared error for a given unkVec"""
function resids(x::Vector{T}) where {T}
    #TODO add weights etc.
    ytrue, tps, expVec, ligVec = getyVec()
    yhat = similar(ytrue, T)

    increasing = tps[2:end] .> tps[1:(end - 1)]

    i = 1
    while i < length(tps)
        jj = 1
        while increasing[i + jj - 1] && i + jj - 1 < length(increasing)
            jj += 1
        end

        # Fix the end case
        if (i + jj - 1) == length(tps) - 1
            jj += 1
        end

        idxs = i:(i + jj - 1)
        vector = vec(fitParams(ligVec[i, 1:3], x, expVec[i, 1:5]))
        yhat[idxs] = runCkinePSTAT(tps[idxs], vector)

        i += jj
    end

    return norm(yhat - ytrue)
end


""" Gets inital unkowns, optimizes them, and returns parameters of best fit"""
function runFit(; itern = 1000000)
    unkVecInit = getUnkVec()

    opts = Optim.Options(outer_iterations = 2, iterations = itern, show_trace = true)
    # TODO: Bounds are artificially tight right now
    fit = optimize(resids, unkVecInit * 0.75, unkVecInit * 1.5, unkVecInit, Fminbox(GradientDescent()), opts, autodiff = :forward)

    return fit.minimizer
end
