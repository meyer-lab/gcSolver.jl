""" This file builds the depletion manuscript, Figure 7. """

using Plots;
using gcSolver;
using DataFrames;
using Statistics;
using GaussianProcesses
plt = Plots;

# Plot of dose response curves
function gpPlot(ligandName, gp, biv = true)

    responseDF = gcSolver.importData()

    time = [0.5, 1, 2, 4]
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Bivalent"] .== biv, filtFrame)

    cellTypeArr = zeros(4)
    cellTypeArr = ["Treg", "CD8", "NK", "Thelper"]

    colors = ["aqua", "coral", "darkorchid", "goldenrod"]
    pl1 = plt.plot()

    μs3d = zeros(length(time), length(doseVec), length(cellTypeArr))
    σ²s3d = similar(μs3d)

    for (i, ITtime) in enumerate(time)
        for (ind, val) in enumerate(cellTypeArr)
            #println("val = ", val)
            filtFrameCell = filter(row -> row["Cell"] .== val, filtFrame)
            fullDataX = filtFrameCell[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc, :Bivalent]]
            
            μs = zeros(length(time), length(doseVec))
            σ²s = similar(μs)

            intrinsLevels = identity.(convert(Matrix, fullDataX)[1, 3:10])
            append!(intrinsLevels, gcSolver.cellHotEnc(val))

            xMat = zeros(length(doseVec), length(intrinsLevels) + 2)
            xMat[:, 1] .= log10.(doseVec)
            xMat[:, 2] .= ITtime
            xMat[:, 3:size(xMat, 2)] .= repeat(intrinsLevels, outer = [1, length(doseVec)])'
            xMat[:, 3] .= log10.(xMat[:, 3])
            xMat[:, 4] .= log10.(xMat[:, 4])
            μs[i, :], σ²s[i, :] = predict_f(gp, xMat')

            #this may need some debugging
            
            for x1 = 1:length(doseVec)
                μs3d[i, x1, ind] = μs[i, x1]
            end
            #println("[i, μs[i, :], ind] = ", [μs[i, :]])
            println("ind μs3d = ", ind, μs3d[:, :, ind])

            for x2 = 1:length(doseVec)
                σ²s3d[i, x2, ind] = σ²s[i, x2]
            end
            #σ²s3d[i, :, ind] = [i, σ²s[i, :], ind]
        end
   end

    μs = zeros(length(time), length(doseVec))
    σ²s = similar(μs3d)

    μs = μs3d[:,:,1] / (μs3d[:,:,2] + μs3d[:,:,3] + μs3d[:,:,4])
    σ²s = σ²s3d[:,:,1] / (σ²s3d[:,:,2] + σ²s3d[:,:,3] + σ²s3d[:,:,4])

    for (i, ITtime) in enumerate(time)
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
                title = string("Treg specificity for Bivalent ", ligandName),
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
            title = string("Treg specificity for Monovalent ", ligandName),
            titlefontsize = 9,
        )
        end
        plt.plot!(doseVec, μs[i, :], c = colors[i], xscale = :log10, ylims = (0,5), yticks = 0:0.5:5, label = ITtime, legend = :bottomright, legendfontsize = 5, markersize = 5)
    
    end

    ylabel!("Treg pSTAT Specificity", yguidefontsize = 7)
    xlabel!("Dose (nM)", xguidefontsize = 7)

    return pl1
end



"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ7()
    l = @layout [a b c d; e f g h; i j k l; m n o p; q r s t; u v w x; z aa bb cc]
    X, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(X', y)
    p1 = gpPlot("IL2", trainedGP, false)
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
