
""" Import all of the data into one dataframe. """
function importData()
    dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")

    yData = DataFrame!(CSV.File(joinpath(dataDir, "WTMuteinsMoments.csv")))
    affDF = DataFrame!(CSV.File(joinpath(dataDir, "mutAffData.csv")))
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
    fitVec = DataFrame!(CSV.File(joinpath(dataDir, "fitTry.csv")))
    return fitVec
end
