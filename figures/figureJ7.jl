""" Treg specificity vs. dose curve (GP model), Figure 7. """

using Plots;
using gcSolver;
using DataFrames;
using Statistics;
using GaussianProcesses
plt = Plots;

# Plot of dose response curves
function specPlot(ligandName, gp, biv = true)

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
    print()

    for (i, ITtime) in enumerate(time)
        for (ind, val) in enumerate(cellTypeArr)
            #println("val = ", val)
            filtFrameCell = filter(row -> row["Cell"] .== val, filtFrame)
            fullDataX = filtFrameCell[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc, :Bivalent]]

            μs = zeros(length(time), length(doseVec))

            intrinsLevels = identity.(convert(Matrix, fullDataX)[1, 3:10])
            append!(intrinsLevels, gcSolver.cellHotEnc(val))

            xMat = zeros(length(doseVec), length(intrinsLevels) + 2)
            xMat[:, 1] .= log10.(doseVec)
            xMat[:, 2] .= ITtime
            xMat[:, 3:size(xMat, 2)] .= repeat(intrinsLevels, outer = [1, length(doseVec)])'
            xMat[:, 3] .= log10.(xMat[:, 3])
            xMat[:, 4] .= log10.(xMat[:, 4])
            μs[i, :] = predict_f(gp, xMat')[1]

            #this may need some debugging

            μs3d[i, :, ind] = 10.0 .^ μs[i, :]

        end
    end

    μs = zeros(length(time), length(doseVec))
    μs = μs3d[:, :, 1] ./ ((μs3d[:, :, 2] .+ μs3d[:, :, 3] .+ μs3d[:, :, 4]) ./ 3)

    for (i, ITtime) in enumerate(time)
        if biv == true
            plt.plot!(
                doseVec,
                μs[i, :],
                c = colors[i],
                xscale = :log10,
                ylims = (0, 5),
                yticks = 0:0.5:5,
                label = "",
                title = string("Treg specificity for Bivalent ", ligandName),
                titlefontsize = 9,
            )
        else
            plt.plot!(
                doseVec,
                μs[i, :],
                c = colors[i],
                xscale = :log10,
                ylims = (0, 5),
                yticks = 0:0.5:5,
                label = "",
                title = string("Treg specificity for Monovalent ", ligandName),
                titlefontsize = 9,
            )
        end
        plt.plot!(
            doseVec,
            μs[i, :],
            c = colors[i],
            xscale = :log10,
            ylims = (0, 30),
            yticks = 0:1:30,
            label = ITtime,
            legend = :bottomright,
            legendfontsize = 5,
            markersize = 5,
        )

    end

    ylabel!("Treg pSTAT Specificity", yguidefontsize = 7)
    xlabel!("Dose (nM)", xguidefontsize = 7)

    return pl1
end



"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ7()
    X, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(X', y)
    p1 = specPlot("IL2", trainedGP, false)
    p2 = specPlot("IL15", trainedGP, false)
    p3 = specPlot("WT N-term", trainedGP, false)
    p4 = specPlot("WT C-term", trainedGP, false)
    p5 = specPlot("V91K C-term", trainedGP, false)
    p6 = specPlot("R38Q N-term", trainedGP, false)
    p7 = specPlot("F42Q N-Term", trainedGP, false)
    p8 = specPlot("N88D C-term", trainedGP, false)
    p9 = specPlot("WT N-term", trainedGP, true)
    p10 = specPlot("R38Q N-term", trainedGP, true)
    p11 = specPlot("H16N N-term", trainedGP, true)
    p12 = specPlot("R38Q/H16N", trainedGP, true)

    ffig = plt.plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, layout = (3, 4), size = (2400, 3200))
    plt.savefig(ffig, joinpath(dirname(pathof(gcSolver)), "..", "figureJ7.svg"))
end
