""" This file builds the depletion manuscript, Figure 7. """

using Plots;
using gcSolver;
using DataFrames;
using Statistics;
using GaussianProcesses
plt = Plots;

# Plot of dose response curves
function gpPlot(ligandName, cellType, gp, biv = true, compType = "none")

    responseDF = gcSolver.importData()

    time = [0.5, 1, 2, 4]
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

    filtFrameTreg = filter(row -> row["Ligand"] .== ligandName, responseDF)


#SHOULD I GET RID OF FILTER FOR CELLTYPE HERE??????


    filter!(row -> row["Cell"] .== "Treg", filtFrameTreg)
    filter!(row -> row["Bivalent"] .== biv, filtFrameTreg)
    
    """filtFrameCD8 = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== "CD8"," filtFrameCD8)
    filter!(row -> row["Bivalent"] .== biv, filtFrameCD8)

    filtFrameNK = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== "NK", filtFrameNK)
    filter!(row -> row["Bivalent"] .== biv, filtFrameNK)

    filtFrameTH = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== "Thelper", filtFrameTH)
    filter!(row -> row["Bivalent"] .== biv, filtFrameTH)"""



    #println("filtFrameTreg = ", filtFrameTreg[1:10, :])

    fullDataX = filtFrameTreg[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc, :Bivalent]]

    #println("fullDataX = ", fullDataX[1:10, :])

    """intrinsLevels = identity.(convert(Matrix, fullDataX)[1, 3:10])
    append!(intrinsLevels, gcSolver.cellHotEnc(cellType))
    xMat = zeros(length(doseVec), length(intrinsLevels) + 2)
    """

    μs = zeros(length(time), length(doseVec))
    σ²s = similar(μs)

    intrinsLevelsR = identity.(convert(Matrix, fullDataX)[1, 3:10])
    append!(intrinsLevelsR, gcSolver.cellHotEnc("Treg"))
    intrinsLevels8 = identity.(convert(Matrix, fullDataX)[1, 3:10])
    append!(intrinsLevels8, gcSolver.cellHotEnc("CD8"))
    intrinsLevelsNK = identity.(convert(Matrix, fullDataX)[1, 3:10])
    append!(intrinsLevelsNK, gcSolver.cellHotEnc("NK"))
    intrinsLevelsH = identity.(convert(Matrix, fullDataX)[1, 3:10])
    append!(intrinsLevelsH, gcSolver.cellHotEnc("Thelper"))
    xMatR = zeros(length(doseVec), length(intrinsLevelsR) + 2)
    xMat8 = zeros(length(doseVec), length(intrinsLevels8) + 2)
    xMatNK = zeros(length(doseVec), length(intrinsLevelsNK) + 2)
    xMatH = zeros(length(doseVec), length(intrinsLevelsH) + 2)

    colors = ["aqua", "coral", "darkorchid", "goldenrod"]
    μsR = zeros(length(time), length(doseVec))
    σ²sR = similar(μsR)
    μs8 = zeros(length(time), length(doseVec))
    σ²s8 = similar(μs8)
    μsNK = zeros(length(time), length(doseVec))
    σ²sNK = similar(μsNK)
    μsH = zeros(length(time), length(doseVec))
    σ²sH = similar(μsH)
    μs3 = zeros(length(time), length(doseVec))
    σ²s3 = similar(μsH)

    pl1 = plt.plot()
 
    for (i, ITtime) in enumerate(time)
        """xMat = zeros(length(doseVec), length(intrinsLevels) + 2)
        xMat[:, 1] .= log10.(doseVec)
        xMat[:, 2] .= ITtime
        xMat[:, 3:size(xMat, 2)] .= repeat(intrinsLevels, outer = [1, length(doseVec)])'
        xMat[:, 3] .= log10.(xMat[:, 3])
        xMat[:, 4] .= log10.(xMat[:, 4])
        μs[i, :], σ²s[i, :] = predict_f(gp, xMat')"""

        xMatR = zeros(length(doseVec), length(intrinsLevelsR) + 2)
        xMatR[:, 1] .= log10.(doseVec)
        xMatR[:, 2] .= ITtime
        xMatR[:, 3:size(xMatR, 2)] .= repeat(intrinsLevelsR, outer = [1, length(doseVec)])'
        xMatR[:, 3] .= log10.(xMatR[:, 3])
        xMatR[:, 4] .= log10.(xMatR[:, 4])
        μsR[i, :], σ²sR[i, :] = predict_f(gp, xMatR')

        xMat8 = zeros(length(doseVec), length(intrinsLevels8) + 2)
        xMat8[:, 1] .= log10.(doseVec)
        xMat8[:, 2] .= ITtime
        xMat8[:, 3:size(xMat8, 2)] .= repeat(intrinsLevels8, outer = [1, length(doseVec)])'
        xMat8[:, 3] .= log10.(xMat8[:, 3])
        xMat8[:, 4] .= log10.(xMat8[:, 4])
        μs8[i, :], σ²s8[i, :] = predict_f(gp, xMat8')

        xMatNK = zeros(length(doseVec), length(intrinsLevelsNK) + 2)
        xMatNK[:, 1] .= log10.(doseVec)
        xMatNK[:, 2] .= ITtime
        xMatNK[:, 3:size(xMatNK, 2)] .= repeat(intrinsLevelsNK, outer = [1, length(doseVec)])'
        xMatNK[:, 3] .= log10.(xMatNK[:, 3])
        xMatNK[:, 4] .= log10.(xMatNK[:, 4])
        μsNK[i, :], σ²sNK[i, :] = predict_f(gp, xMatNK')

        xMatH = zeros(length(doseVec), length(intrinsLevelsH) + 2)
        xMatH[:, 1] .= log10.(doseVec)
        xMatH[:, 2] .= ITtime
        xMatH[:, 3:size(xMatH, 2)] .= repeat(intrinsLevelsH, outer = [1, length(doseVec)])'
        xMatH[:, 3] .= log10.(xMatH[:, 3])
        xMatH[:, 4] .= log10.(xMatH[:, 4])
        μsH[i, :], σ²sH[i, :] = predict_f(gp, xMatH')

        μs3[i, :] = μs8[i, :] + μsNK[i, :] + μsH[i, :]
        println(μs3[1, 1], " = ", μs8[1, 1], " + ", μsNK[1, 1], " + ", μsH[1, 1])
        σ²s3[i, :] = σ²s8[i, :] + σ²sNK[i, :] + σ²sH[i, :]

        println(μs[i, :], " = ", μsR[i, :], " / ", μs3[1, :])
        μs[i, :] = μsR[i, :] / μs3[i, :]
        σ²s[i, :] = σ²sR[i, :] / σ²s3[i, :]

        #println("y = ", μs[i, :])
        #println("σ²s = ", σ²s[1:4, :])

        if biv == true
            plt.plot!(
                doseVec,
                [μs[i, :] μs[i, :]],
                fillrange = [μs[i, :] .- σ²s[i, :] μs[i, :] .+ σ²s[i, :]],
                fillalpha = 0.3,
                c = colors[i],
                xscale = :log10,
                ylims = (0,5),
                yticks = 0:0.5:5,
                label = "",
                title = string(cellType, " Response to Bivalent ", ligandName, " GP Model"),
                titlefontsize = 9,
            )
        else
            plt.plot!(
            doseVec,
            [μs[i, :] μs[i, :]],
            fillrange = [μs[i, :] .- σ²s[i, :] μs[i, :] .+ σ²s[i, :]],
            fillalpha = 0.3,
            c = colors[i],
            xscale = :log10,
            ylims = (0,5),
            yticks = 0:0.5:5,
            label = "",
            title = string(cellType, " Response to Monovalent ", ligandName, " GP Model"),
            titlefontsize = 9,
        )
        end
        plt.plot!(doseVec, μs[i, :], c = colors[i], xscale = :log10, ylims = (0,5), yticks = 0:0.5:5, label = ITtime, legend = :bottomright, legendfontsize = 5, markersize = 5)

    end

    if compType != "none"
        intrinsLevelsComp = identity.(convert(Matrix, fullDataX)[1, 3:10])
        append!(intrinsLevelsComp, gcSolver.cellHotEnc(compType))
        xMatComp = zeros(length(doseVec), length(intrinsLevelsComp) + 2)

        colorsComp = ["turquoise3", "coral3", "darkorchid4", "darkgoldenrod"]
        μsComp = zeros(length(time), length(doseVec))
        σ²sComp = similar(μs)

        #pl1 = plt.plot()

        for (i, ITtime) in enumerate(time)
            xMatComp = zeros(length(doseVec), length(intrinsLevelsComp) + 2)
            xMatComp[:, 1] .= log10.(doseVec)
            xMatComp[:, 2] .= ITtime
            xMatComp[:, 3:size(xMatComp, 2)] .= repeat(intrinsLevelsComp, outer = [1, length(doseVec)])'
            xMatComp[:, 3] .= log10.(xMatComp[:, 3])
            xMatComp[:, 4] .= log10.(xMatComp[:, 4])
            μsComp[i, :], σ²sComp[i, :] = predict_f(gp, xMatComp')

            #assuming this adds to current plot
            #linestyle not working, maybe dots are too big and looks like solid line?
            if biv == true
                plt.plot!(
                    doseVec,
                    [μsComp[i, :] μsComp[i, :]],
                    #fillrange = [μsComp[i, :] .- σ²sComp[i, :] μsComp[i, :] .+ σ²sComp[i, :]],
                    #fillalpha = 0.3,
                    c = colorsComp[i],
                    #linestyle = :dot,
                    xscale = :log10,
                    ylims = (0,5),
                    yticks = 0:0.5:5,
                    label = "",
                    title = string(cellType, " Response to Bivalent ", ligandName, " GP Model"),
                    titlefontsize = 9,
                )
            else
                plt.plot!(
                    doseVec,
                    [μsComp[i, :] μsComp[i, :]],
                    #fillrange = [μsComp[i, :] .- σ²sComp[i, :] μsComp[i, :] .+ σ²sComp[i, :]],
                    #fillalpha = 0.3,
                    c = colorsComp[i],
                    #linestyle = :dot,
                    xscale = :log10,
                    ylims = (0,5),
                    yticks = 0:0.5:5,
                    label = "",
                    title = string(cellType, " Response to Monovalent ", ligandName, " GP Model"),
                    titlefontsize = 9,
                )
            end
            labelString = string(compType, " ", ITtime)
            plt.plot!(
                doseVec,
                μsComp[i, :],
                c = colorsComp[i],
                xscale = :log10,
                ylims = (0,5),
                yticks = 0:0.5:5,
                label = labelString,
                legend = :bottomright,
                legendfontsize = 5,
                markersize = 5,
            )
        end
    end

    ylabel!("pSTAT", yguidefontsize = 7)
    xlabel!("Dose (nM)", xguidefontsize = 7)

    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ7()
    l = @layout [a b c d; e f g h; i j k l; m n o p; q r s t; u v w x; z aa bb cc]
    X, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(X', y)
    p1 = gpPlot("IL2", "Treg", trainedGP, false)
    """p1 = gpPlot("IL2", "Treg", trainedGP, false)
    p2 = gpPlot("IL2", "Thelper", trainedGP, false)
    p3 = gpPlot("IL2", "NK", trainedGP, false)
    p4 = gpPlot("IL2", "CD8", trainedGP, false)
    p5 = gpPlot("IL15", "Treg", trainedGP, false)
    p6 = gpPlot("IL15", "Thelper", trainedGP, false)
    p7 = gpPlot("IL15", "NK", trainedGP, false)
    p8 = gpPlot("IL15", "CD8", trainedGP, false)
    p9 = gpPlot("WT N-term", "Treg", trainedGP)
    p10 = gpPlot("WT N-term", "Thelper", trainedGP)
    p11 = gpPlot("WT N-term", "NK", trainedGP)
    p12 = gpPlot("WT N-term", "CD8", trainedGP)
    p13 = gpPlot("WT N-term", "Treg", trainedGP, false)
    p14 = gpPlot("WT N-term", "Thelper", trainedGP, false)
    p15 = gpPlot("WT N-term", "NK", trainedGP, false)
    p16 = gpPlot("WT N-term", "CD8", trainedGP, false)
    p17 = gpPlot("R38Q N-term", "Treg", trainedGP)
    p18 = gpPlot("R38Q N-term", "Thelper", trainedGP)
    p19 = gpPlot("R38Q N-term", "NK", trainedGP)
    p20 = gpPlot("R38Q N-term", "CD8", trainedGP)
    p21 = gpPlot("R38Q N-term", "Treg", trainedGP, false)
    p22 = gpPlot("R38Q N-term", "Thelper", trainedGP, false)
    p23 = gpPlot("R38Q N-term", "NK", trainedGP, false)
    p24 = gpPlot("R38Q N-term", "CD8", trainedGP, false)
    p25 = gpPlot("H16N N-term", "Treg", trainedGP)
    p26 = gpPlot("H16N N-term", "Thelper", trainedGP)
    p27 = gpPlot("H16N N-term", "NK", trainedGP)
    p28 = gpPlot("H16N N-term", "CD8", trainedGP)
    p29 = gpPlot("R38Q/H16N", "Treg", trainedGP)
    p30 = gpPlot("R38Q/H16N", "Thelper", trainedGP)
    p31 = gpPlot("R38Q/H16N", "NK", trainedGP)
    p32 = gpPlot("R38Q/H16N", "CD8", trainedGP)
    p33 = gpPlot("WT C-term", "Treg", trainedGP, false)
    p34 = gpPlot("WT C-term", "Thelper", trainedGP, false)
    p35 = gpPlot("WT C-term", "NK", trainedGP, false)
    p36 = gpPlot("WT C-term", "CD8", trainedGP, false)
    p37 = gpPlot("V91K C-term", "Treg", trainedGP, false)
    p38 = gpPlot("V91K C-term", "Thelper", trainedGP, false)
    p39 = gpPlot("V91K C-term", "NK", trainedGP, false)
    p40 = gpPlot("V91K C-term", "CD8", trainedGP, false)
    p41 = gpPlot("F42Q N-Term", "Treg", trainedGP, false)
    p42 = gpPlot("F42Q N-Term", "Thelper", trainedGP, false)
    p43 = gpPlot("F42Q N-Term", "NK", trainedGP, false)
    p44 = gpPlot("F42Q N-Term", "CD8", trainedGP, false)
    p45 = gpPlot("N88D C-term", "Treg", trainedGP, false)
    p46 = gpPlot("N88D C-term", "Thelper", trainedGP, false)
    p47 = gpPlot("N88D C-term", "NK", trainedGP, false)
    p48 = gpPlot("N88D C-term", "CD8", trainedGP, false)"""
    

    #draw(SVG("figureJ1.svg", 1000px, 800px), p1)
    """ffig = plt.plot(
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
        p25,
        p26,
        p27,
        p28,
        p29,
        p30,
        p31,
        p32,
        p33,
        p34,
        p35,
        p36,
        p37,
        p38,
        p39,
        p40,
        p41,
        p42,
        p43,
        p44,
        p45,
        p46,
        p47,
        p48,        
        layout = (12,4),
        size = (2000, 5000),
    )"""
    ffig = plt.plot(p1, size = (1600, 2400))
    plt.savefig(ffig, joinpath(dirname(pathof(gcSolver)), "..", "figureJ7.svg"))
end
