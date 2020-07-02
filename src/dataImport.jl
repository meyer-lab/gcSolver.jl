
""" Import all of the data into one dataframe. """
@memoize function importData()
	dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")

	yData = DataFrame!(CSV.File(joinpath(dataDir, "WTMuteinsMoments.csv")))
    affDF = DataFrame!(CSV.File(joinpath(dataDir, "mutAffData.csv")))
    exprDF = DataFrame!(CSV.File(joinpath(dataDir, "RecQuantitation.csv")))

    fullData = innerjoin(yData, affDF, on = :Ligand => :Mutein)

    return fullData, exprDF
end

""" Import the saved fit. """
function importFit()
	dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")
	fitVec = CSV.read(joinpath(dataDir, "fitTry.csv"))
	return fitVec
end