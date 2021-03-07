
""" Import a file. """
function importFile(name::String)
    dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")
    return DataFrame(CSV.File(joinpath(dataDir, name)))
end

""" Import all of the data into one dataframe. """
function importData(monomeric = false)
    yData = importFile("WTDimericMutSingleCellData.csv")
    monDF = importFile("MonomericMutSingleCellData.csv")

    append!(yData, monDF)
    if monomeric
        yData = filter(row -> row.Bivalent ∈ [0], yData)
    end

    affDF = importFile("WTmutAffData.csv")
    exprDF = importFile("RecQuantitation.csv")

    exprDF = stack(exprDF, [:Treg, :Thelper, :NK, :CD8]; variable_name = "Cell")
    exprDF = unstack(exprDF, :Receptor, :value)
    select!(exprDF, Not(:Moment))

    fullData = innerjoin(yData, affDF, on = :Ligand => :Mutein)
    fullData = leftjoin(fullData, exprDF, on = :Cell)

    return fullData
end

""" Import the saved fit. """
function importFit()
    return importFile("actualFit.csv")
end

""" Import Date pSTAT conversion factors. """
function importConvFrame()
    return importFile("DateConvFrame.csv")
end

"""Creates Sigma for Var Propagation"""
function getSigma(cellType)
    sigma = zeros(3, 3)

    momentDF = importFile("receptor_moments.csv")
    covDF = importFile("receptor_covariances.csv")

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
