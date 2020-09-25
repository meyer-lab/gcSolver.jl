""" This file builds the depletion manuscript, Figure 5. """

using Plots;
using gcSolver;
using DataFrames;
using Statistics;
plt = Plots;

# Plot of dose response curves
function gpPlotVar(ligandName, cellType, gp, cov = false, biv = true)
    nPoints = 100
    responseDF = gcSolver.importData()
    sigma = gcSolver.getSigma(cellType)
    time = [0.5, 1, 2, 4]
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    filter!(row -> row["Bivalent"] .== biv, filtFrame)

    if cov
        realDataDF = filtFrame[!, [:Dose, :Time, :alphStatCov]]
        realDataDF = groupby(realDataDF, [:Time, :Dose])
        realDataDF = combine(realDataDF, :alphStatCov => mean)
    else
        realDataDF = filtFrame[!, [:Dose, :Time, :Variance]]
        realDataDF = groupby(realDataDF, [:Time, :Dose])
        realDataDF = combine(realDataDF, :Variance => mean)
    end

    fullDataX = filtFrame[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc, :Bivalent]]
    intrinsLevels = identity.(convert(Matrix, fullDataX)[1, 3:10])
    append!(intrinsLevels, gcSolver.cellHotEnc(cellType))
    xMat = zeros(length(doseVec), length(intrinsLevels) + 2)

    colors = ["aqua", "coral", "darkorchid", "goldenrod"]
    vars = zeros(nPoints)

    pl1 = plt.plot()
    doseVecPred = exp10.(range(log10.(doseVec[1]), log10.(doseVec[length(doseVec)]), length = nPoints))

    for (i, ITtime) in enumerate(time)
        xMat = zeros(nPoints, length(intrinsLevels) + 2)
        xMat[:, 1] .= log10.(doseVecPred)
        xMat[:, 2] .= ITtime
        xMat[:, 3:size(xMat, 2)] .= repeat(intrinsLevels, outer = [1, nPoints])'
        xMat[:, 3] .= log10.(xMat[:, 3])
        xMat[:, 4] .= log10.(xMat[:, 4])
        for ii in (1:nPoints)
            vars[ii] = log10.(gcSolver.runCkineVarPropGP(gp, xMat[ii, :], sigma, cov))[1]
        end

        if cov
            plt.plot!(doseVecPred, vars, c = colors[i], xscale = :log10, label = ITtime, legend = :bottomright, legendfontsize = 5, markersize = 5)
            if length(log10.(abs.(realDataDF[realDataDF.Time .== ITtime, :].alphStatCov_mean) .+ 1)) > 0
                absCovs = abs.(realDataDF[realDataDF.Time .== ITtime, :].alphStatCov_mean)
                plt.scatter!(
                    doseVec,
                    log10.(absCovs) .+ 1,
                    c = colors[i],
                    xscale = :log10,
                    label = "",
                    title = string(cellType, " Response to ", ligandName, " GP Model"),
                    titlefontsize = 9,
                )
            end

        else
            plt.plot!(doseVecPred, vars, c = colors[i], xscale = :log10, label = ITtime, legend = :bottomright, legendfontsize = 5, markersize = 5)
            if length(log10.(realDataDF[realDataDF.Time .== ITtime, :].Variance_mean .+ 1)) > 0
                plt.scatter!(
                    doseVec,
                    log10.(realDataDF[realDataDF.Time .== ITtime, :].Variance_mean) .+ 1,
                    c = colors[i],
                    xscale = :log10,
                    label = "",
                    title = string(cellType, " Response to ", ligandName, " GP Model"),
                    titlefontsize = 9,
                )
            end
        end

    end

    if cov
        ylabel!("log pSTAT/IL2Ra Covariance", yguidefontsize = 7)
        xlabel!("Dose (nM)", xguidefontsize = 7)
    else
        ylabel!("log pSTAT Variance", yguidefontsize = 7)
        xlabel!("Dose (nM)", xguidefontsize = 7)
    end

    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ5()
    l = @layout [a b c d; e f g h; i j k l; m n o p; q r s t; u v w x; z aa]
    X, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(X', y)
    p1 = gpPlotVar("IL2", "Treg", trainedGP, true, false)
    p2 = gpPlotVar("IL2", "Thelper", trainedGP, true, false)
    p3 = gpPlotVar("IL2", "NK", trainedGP, false, false)
    p4 = gpPlotVar("IL2", "CD8", trainedGP, false, false)
    p5 = gpPlotVar("IL15", "Treg", trainedGP, true, false)
    p6 = gpPlotVar("IL15", "Thelper", trainedGP, true, false)
    p7 = gpPlotVar("IL15", "NK", trainedGP, false, false)
    p8 = gpPlotVar("IL15", "CD8", trainedGP, false, false)
    p9 = gpPlotVar("R38Q/H16N", "Treg", trainedGP, true)
    p10 = gpPlotVar("R38Q/H16N", "Thelper", trainedGP, true)
    p11 = gpPlotVar("R38Q/H16N", "NK", trainedGP)
    p12 = gpPlotVar("R38Q/H16N", "CD8", trainedGP)
    p13 = gpPlotVar("WT N-term", "Treg", trainedGP, true)
    p14 = gpPlotVar("WT N-term", "Thelper", trainedGP, true)
    p15 = gpPlotVar("WT N-term", "NK", trainedGP)
    p16 = gpPlotVar("WT N-term", "CD8", trainedGP)
    p17 = gpPlotVar("H16N N-term", "Treg", trainedGP, true)
    p18 = gpPlotVar("H16N N-term", "Thelper", trainedGP, true)
    p19 = gpPlotVar("H16N N-term", "NK", trainedGP)
    p20 = gpPlotVar("H16N N-term", "CD8", trainedGP)
    p21 = gpPlotVar("R38Q N-term", "Treg", trainedGP, true)
    p22 = gpPlotVar("R38Q N-term", "Thelper", trainedGP, true)
    p23 = gpPlotVar("R38Q N-term", "NK", trainedGP)
    p24 = gpPlotVar("R38Q N-term", "CD8", trainedGP)
    p26 = gpPlotVar("V91K C-term", "CD8", trainedGP, true, false)
    p28 = gpPlotVar("V91K C-term", "CD8", trainedGP, false, false)
    #draw(SVG("figureJ1.svg", 1000px, 800px), p1)
    ffig = plt.plot(
        p1,
        p2,
        p3,
        p4,
        p5,
        p6,
        p7,
        p8,
        p9,
        p10,
        p11,
        p12,
        p13,
        p14,
        p15,
        p16,
        p17,
        p18,
        p19,
        p20,
        p21,
        p22,
        p23,
        p24,
        p26,
        p28,
        layout = l,
        size = (1600, 2400),
    )
    plt.savefig(ffig, joinpath(dirname(pathof(gcSolver)), "..", "figureJ5.svg"))
end
