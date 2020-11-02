""" This file builds the depletion manuscript, Figure 8. """

using Plots;
using gcSolver;
using DataFrames;
using Statistics;
using GaussianProcesses
plt = Plots;

# Plot of dose response curves
function gpPlot(ligandName, cellType, gp, time)

    responseDF = gcSolver.importData()

    valency = [1, 0]
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    filter!(row -> row["Time"] .== time, filtFrame)
    #filter!(row -> string(row["Date"]) .== date, filtFrame)

    colors = ["aqua", "goldenrod"]

    pl1 = plt.plot()

    for (ind, val) in enumerate(valency)
        valDF = filter(row -> row["Bivalent"] .== val, filtFrame)

        valDataDF = valDF[!, [:Dose, :Bivalent, :Mean]]
        valDataDF = groupby(valDataDF, [:Bivalent, :Dose])
        valDataDF = combine(valDataDF, :Mean => mean)

        fullDataX = valDF[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc, :Bivalent]]

        intrinsLevels = identity.(convert(Matrix, fullDataX)[1, 3:10])
        append!(intrinsLevels, gcSolver.cellHotEnc(cellType))
        xMat = zeros(length(doseVec), length(intrinsLevels) + 2)

        μs = zeros(length(valency), length(doseVec))
        #println("valLength = ", length(valency))
        σ²s = similar(μs)

        pl1 = plt.plot()
        xMat = zeros(length(doseVec), length(intrinsLevels) + 2)
        xMat[:, 1] .= log10.(doseVec)
        xMat[:, 2] .= time
        xMat[:, 3:size(xMat, 2)] .= repeat(intrinsLevels, outer = [1, length(doseVec)])'
        xMat[:, 3] .= log10.(xMat[:, 3])
        xMat[:, 4] .= log10.(xMat[:, 4])
        μs[ind, :], σ²s[ind, :] = predict_f(gp, xMat')

        plt.plot!(
        doseVec,
        [μs[ind, :] μs[ind, :]],
        fillrange = [μs[ind, :] .- σ²s[ind, :] μs[ind, :] .+ σ²s[ind, :]],
        fillalpha = 0.3,
        c = colors[ind],
        xscale = :log10,
        ylims = (0,5),
        yticks = 0:0.5:5,
        label = "",
        title = string(cellType, " Response to ", ligandName, " at t = ", time),
        titlefontsize = 9)

        plt.plot!(doseVec, μs[ind, :], c = colors[ind], xscale = :log10, ylims = (0,5), yticks = 0:0.5:5, label = val, legend = :bottomright, legendfontsize = 5, markersize = 5)

        if length(log10.(valDataDF[valDataDF.Bivalent .== val, :].Mean_mean .+ 1)) > 0
            plt.scatter!(doseVec, log10.(valDataDF[valDataDF.Bivalent .== val, :].Mean_mean .+ 1), c = colors[ind], xscale = :log10, label = "")
        end
    end

    ylabel!("pSTAT", yguidefontsize = 7)
    xlabel!("Dose (nM)", xguidefontsize = 7)

    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ8()
    l = @layout [a b c d; e f g h; i j k l; m n o p; q r s t; u v w x; z aa bb cc]
    X, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(X', y)
    p1 = gpPlot("R38Q N-term", "Treg", trainedGP, 2)
    p2 = gpPlot("WT N-term", "Treg", trainedGP, 2)
    
    ffig = plt.plot(p1, p2, layout = (2,1), size = (1600, 2400))
    plt.savefig(ffig, joinpath(dirname(pathof(gcSolver)), "..", "figureJ8.svg"))
end
