""" This file builds the depletion manuscript, Figure 1. """

using Plots;
plt = Plots;

# Plot of dose response curves
function gpPlot(ligandName, cellType, gp, compType = "none")
    responseDF = importData()
    time = [0.5, 1, 2, 4]
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    #filter!(row -> string(row["Date"]) .== date, filtFrame)

    realDataDF = filtFrame[!, [:Dose, :Time, :Mean]]
    realDataDF = groupby(realDataDF, [:Time, :Dose])
    realDataDF = combine(realDataDF, :Mean => mean)

    fullDataX = filtFrame[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc]]
    
    intrinsLevels = identity.(convert(Matrix, fullDataX)[1, 3:9])
    append!(intrinsLevels, cellHotEnc(cellType))
    xMat = zeros(length(doseVec), length(intrinsLevels) + 2)

    colors = ["aqua", "coral", "darkorchid", "goldenrod"]
    μs = zeros(length(time), length(doseVec))
    σ²s = similar(μs)

    pl1 = plt.plot()

    for (i, ITtime) in enumerate(time)
        xMat = zeros(length(doseVec), length(intrinsLevels) + 2)
        xMat[:, 1] .= log10.(doseVec)
        xMat[:, 2] .= ITtime
        xMat[:, 3:size(xMat,2)] .= repeat(intrinsLevels, outer = [1, length(doseVec)])'
        xMat[:, 3] .= log10.(xMat[:, 3])
        xMat[:, 4] .= log10.(xMat[:, 4])
        μs[i, :], σ²s[i, :] = predict_f(gp, xMat')
        plt.plot!(
            doseVec,
            [μs[i, :] μs[i, :]],
            fillrange = [μs[i, :] .- σ²s[i, :] μs[i, :] .+ σ²s[i, :]],
            fillalpha = 0.3,
            c = colors[i],
            xscale = :log10,
            label = "",
            title = string(cellType, " Response to ", ligandName, " GP Model"),
            titlefontsize = 9,
        )
        plt.plot!(doseVec, μs[i, :], c = colors[i], xscale = :log10, label = ITtime, legend = :bottomright, legendfontsize = 5, markersize = 5)

        if length(log10.(realDataDF[realDataDF.Time .== ITtime, :].Mean_mean .+ 1)) > 0
            plt.scatter!(doseVec, log10.(realDataDF[realDataDF.Time .== ITtime, :].Mean_mean .+ 1), c = colors[i], xscale = :log10, label = "")
        end

    end
    
    if compType != "none"
        intrinsLevelsComp = identity.(convert(Matrix, fullDataX)[1, 3:9])
        append!(intrinsLevelsComp, cellHotEnc(compType))
        xMatComp = zeros(length(doseVec), length(intrinsLevelsComp) + 2)

        colorsComp = ["turquoise3", "coral3", "darkorchid4", "darkgoldenrod"]
        μsComp = zeros(length(time), length(doseVec))
        σ²sComp = similar(μs)

        #pl1 = plt.plot()

        for (i, ITtime) in enumerate(time)
            xMatComp = zeros(length(doseVec), length(intrinsLevelsComp) + 2)
            xMatComp[:, 1] .= log10.(doseVec)
            xMatComp[:, 2] .= ITtime
            xMatComp[:, 3:size(xMatComp,2)] .= repeat(intrinsLevelsComp, outer = [1, length(doseVec)])'
            xMatComp[:, 3] .= log10.(xMatComp[:, 3])
            xMatComp[:, 4] .= log10.(xMatComp[:, 4])
            μsComp[i, :], σ²sComp[i, :] = predict_f(gp, xMatComp')
            
            #assuming this adds to current plot
            #linestyle not working, maybe dots are too big and looks like solid line?
            plt.plot!(
                doseVec,
                [μsComp[i, :] μsComp[i, :]],
                #fillrange = [μsComp[i, :] .- σ²sComp[i, :] μsComp[i, :] .+ σ²sComp[i, :]],
                #fillalpha = 0.3,
                c = colorsComp[i],
                #linestyle = :dot,
                xscale = :log10,
                label = "",
                title = string(cellType, " Response to ", ligandName, " GP Model"),
                titlefontsize = 9,
            )
            labelString = string(compType, " ", ITtime)
            plt.plot!(doseVec, μsComp[i, :], c = colorsComp[i], xscale = :log10, label = labelString, legend = :bottomright, legendfontsize = 5, markersize = 5)

            # checks that we have data before attempting to plot it
            if length(log10.(realDataDF[realDataDF.Time .== ITtime, :].Mean_mean .+ 1)) > 0
                plt.scatter!(doseVec, log10.(realDataDF[realDataDF.Time .== ITtime, :].Mean_mean .+ 1), c = colorsComp[i], xscale = :log10, label = "")
            end
        end
    end

    ylabel!("pSTAT", yguidefontsize = 7)
    xlabel!("Dose (nM)", xguidefontsize = 7)

    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ3()
    l = @layout [a b c d; e f g h; i j k l; m n o p; q r s t; u v w x]
    X, y, df = getGPdata()
    trainedGP = gaussianProcess(X', y)
    #p1 = gpPlot("IL2", "Treg", trainedGP)
    p1 = gpPlot("IL2", "Treg", trainedGP, "NK")
    p2 = gpPlot("IL2", "Thelper", trainedGP)
    p3 = gpPlot("IL2", "NK", trainedGP)
    p4 = gpPlot("IL2", "CD8", trainedGP)
    p5 = gpPlot("IL15", "Treg", trainedGP)
    p6 = gpPlot("IL15", "Thelper", trainedGP)
    p7 = gpPlot("IL15", "NK", trainedGP)
    p8 = gpPlot("IL15", "CD8", trainedGP)
    p9 = gpPlot("R38Q/H16N", "Treg", trainedGP)
    p10 = gpPlot("R38Q/H16N", "Thelper", trainedGP)
    p11 = gpPlot("R38Q/H16N", "NK", trainedGP)
    p12 = gpPlot("R38Q/H16N", "CD8", trainedGP)
    p13 = gpPlot("WT N-term", "Treg", trainedGP)
    p14 = gpPlot("WT N-term", "Thelper", trainedGP)
    p15 = gpPlot("WT N-term", "NK", trainedGP)
    p16 = gpPlot("WT N-term", "CD8", trainedGP)
    p17 = gpPlot("H16N N-term", "Treg", trainedGP)
    p18 = gpPlot("H16N N-term", "Thelper", trainedGP)
    p19 = gpPlot("H16N N-term", "NK", trainedGP)
    p20 = gpPlot("H16N N-term", "CD8", trainedGP)
    p21 = gpPlot("R38Q N-term", "Treg", trainedGP)
    p22 = gpPlot("R38Q N-term", "Thelper", trainedGP)
    p23 = gpPlot("R38Q N-term", "NK", trainedGP)
    p24 = gpPlot("R38Q N-term", "CD8", trainedGP)
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
        layout = l,
        size = (1600, 2400),
    )
    #ffig = plt.plot(p1, size = (1600, 2400))
    plt.savefig(ffig, joinpath(dirname(pathof(gcSolver)), "..", "figureJ3.svg"))
end
