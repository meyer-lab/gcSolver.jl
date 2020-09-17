
""" Import all of the data into one dataframe. """
function importData(monomeric = false)
    dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")

    yData = DataFrame!(CSV.File(joinpath(dataDir, "WTMuteinsMoments.csv")))
    if monomeric
        yData = filter(row -> row.Ligand ∈ ["IL2", "IL15"], yData)
        monDF = DataFrame!(CSV.File(joinpath(dataDir, "MonomericMutpSTAT.csv")))
        append!(yData, monDF)
    end
    affDF = DataFrame!(CSV.File(joinpath(dataDir, "WTmutAffData.csv")))
    exprDF = DataFrame!(CSV.File(joinpath(dataDir, "RecQuantitation.csv")))

    exprDF = stack(exprDF, [:Treg, :Thelper, :NK, :CD8]; variable_name = "Cell")
    exprDF = unstack(exprDF, :Receptor, :value)
    select!(exprDF, Not(:Moment))

    fullData = innerjoin(yData, affDF, on = :Ligand => :Mutein)
    fullData = leftjoin(fullData, exprDF, on = :Cell)

    return fullData
end

""" Import the saved fit. """
function importFit()
    dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")
    fitVec = DataFrame!(CSV.File(joinpath(dataDir, "MonomericFit.csv")))
    return fitVec
end

""" Import Date pSTAT conversion factors. """
function importConvFrame()
    dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")
    convFrame = DataFrame!(CSV.File(joinpath(dataDir, "DateConvFrame.csv")))
    return convFrame
end

"""Creates Sigma for Var Propagation"""
function getSigma(cellType)
    dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")
    sigma = zeros(3, 3)

    momentDF = DataFrame!(CSV.File(joinpath(dataDir, "receptor_moments.csv")))
    covDF = DataFrame!(CSV.File(joinpath(dataDir, "receptor_covariances.csv")))

    momentDF = momentDF[!, ["Cell Type", "Receptor", "Variance"]]
    momentDF = groupby(momentDF, ["Cell Type", "Receptor"])
    momentDF = combine(momentDF, :Variance => mean)
    filter!(row -> row["Cell Type"] .== cellType, momentDF)
    for i = 1:3
        sigma[i, i] = momentDF[i, "Variance_mean"]
    end

    covDF = covDF[!, ["Cell Type", "CD25:Receptor", "Covariance"]]
    covDF = groupby(covDF, ["Cell Type", "CD25:Receptor"])
    covDF = combine(covDF, :Covariance => mean)
    filter!(row -> row["Cell Type"] .== "Treg", covDF)
    for i = 2:3
        sigma[1, i] = sigma[i, 1] = covDF[i - 1, "Covariance_mean"]
    end

    return sigma
end
