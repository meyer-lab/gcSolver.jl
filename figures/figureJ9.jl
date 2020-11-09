""" This file builds the depletion manuscript, Figure 9. """

using Plots;
using gcSolver;
using DataFrames;
using Statistics;
using GaussianProcesses
plt = Plots;

# Plot of dose response curves
function gpPlot(ligandName, cellType, gp, time)

    responseDF = gcSolver.importData()

    valency = [true, false]
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    filter!(row -> row["Time"] .== time, filtFrame)
    #filter!(row -> string(row["Date"]) .== date, filtFrame)

    colors = ["aqua", "goldenrod"]
    leg = ["Bivalent", "Monovalent"]

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

        μs = zeros(length(doseVec))
        #println("valLength = ", length(valency))
        σ²s = similar(μs)

        xMat = zeros(length(doseVec), length(intrinsLevels) + 2)
        xMat[:, 1] .= doseVec
        xMat[:, 2] .= time
        xMat[:, 3:size(xMat, 2)] .= repeat(intrinsLevels, outer = [1, length(doseVec)])'
        xMat[:, 3] .= xMat[:, 3]
        xMat[:, 4] .= xMat[:, 4]
        μs[:], σ²s[:] = predict_f(gp, xMat')

        """plt.plot!(
        doseVec,
        [μs μs],
        fillrange = [μs .- σ²s μs .+ σ²s],
        fillalpha = 0.3,
        c = colors[ind],
        #ylims = (0,5),
        #yticks = 0:0.5:5,
        label = "",
        title = string(cellType, " Response to ", ligandName, " at t = ", time),
        titlefontsize = 9)"""

        plt.plot!(doseVec, μs, c = colors[ind], 
        #ylims = (0,5), 
        #yticks = 0:0.5:5, 
        label = leg[ind], title = string(cellType, " Response to ", ligandName, " at t = ", time), legend = :bottomright, legendfontsize = 5, markersize = 5)

        if length(valDataDF[valDataDF.Bivalent .== val, :].Mean_mean .+ 1) > 0
            plt.scatter!(doseVec, valDataDF[valDataDF.Bivalent .== val, :].Mean_mean .+ 1, c = colors[ind], label = "")
        end
    end

    ylabel!("pSTAT", yguidefontsize = 7)
    xlabel!("Dose (nM)", xguidefontsize = 7)

    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ9()
    l = @layout [a b c d; e f g h; i j k l; m n o p; q r s t; u v w x; z aa bb cc]
    X, y, df = gcSolver.getGPdataNL()
    trainedGP = gcSolver.gaussianProcess(X', y)
    p1 = gpPlot("R38Q N-term", "Treg", trainedGP, 2)
    p2 = gpPlot("WT N-term", "Treg", trainedGP, 2)
    p3 = gpPlot("R38Q N-term", "CD8", trainedGP, 2)
    #p4 = gpPlot("WT N-term", "CD8", trainedGP, 2)
    p5 = gpPlot("R38Q N-term", "NK", trainedGP, 2)
    #p6 = gpPlot("WT N-term", "NK", trainedGP, 2)
    p7 = gpPlot("R38Q N-term", "Thelper", trainedGP, 2)
    p8 = gpPlot("WT N-term", "Thelper", trainedGP, 2)

    p9 = gpPlot("R38Q N-term", "Treg", trainedGP, 4)
    p10 = gpPlot("WT N-term", "Treg", trainedGP, 4)
    p11 = gpPlot("R38Q N-term", "CD8", trainedGP, 4)
    p12 = gpPlot("WT N-term", "CD8", trainedGP, 4)
    p13 = gpPlot("R38Q N-term", "NK", trainedGP, 4)
    p14 = gpPlot("WT N-term", "NK", trainedGP, 4)
    p15 = gpPlot("R38Q N-term", "Thelper", trainedGP, 4)
    p16 = gpPlot("WT N-term", "Thelper", trainedGP, 4)

    #p17 = gpPlot("R38Q N-term", "Treg", trainedGP, 0.5)
    p18 = gpPlot("WT N-term", "Treg", trainedGP, 0.5)
    #p19 = gpPlot("R38Q N-term", "CD8", trainedGP, 0.5)
    p20 = gpPlot("WT N-term", "CD8", trainedGP, 0.5)
    #p21 = gpPlot("R38Q N-term", "NK", trainedGP, 0.5)
    p22 = gpPlot("WT N-term", "NK", trainedGP, 0.5)
    #p23 = gpPlot("R38Q N-term", "Thelper", trainedGP, 0.5)
    p24 = gpPlot("WT N-term", "Thelper", trainedGP, 0.5)

    p25 = gpPlot("R38Q N-term", "Treg", trainedGP, 1)
    p26 = gpPlot("WT N-term", "Treg", trainedGP, 1)
    p27 = gpPlot("R38Q N-term", "CD8", trainedGP, 1)
    p28 = gpPlot("WT N-term", "CD8", trainedGP, 1)
    p29 = gpPlot("R38Q N-term", "NK", trainedGP, 1)
    p30 = gpPlot("WT N-term", "NK", trainedGP, 1)
    p31 = gpPlot("R38Q N-term", "Thelper", trainedGP, 1)
    p32 = gpPlot("WT N-term", "Thelper", trainedGP, 1)
    
    ffig = plt.plot(p1,
    p2,
    p3,
    #p4,
    p5,
    #p6,
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
    #p17,
    p18,
    #p19,
    p20,
    #p21,
    p22,
    #p23,
    p24,
    p25,
    p26,
    p27,
    p28,
    p29,
    p30,
    p31,
    p32,
    layout = (13,2), size = (1200, 6400))
    #ffig = plt.plot(p1, size = (600,600))
    plt.savefig(ffig, joinpath(dirname(pathof(gcSolver)), "..", "figureJ9.svg"))
    

end
