""" This file builds the depletion manuscript, Figure 1. """

using Plots;
plt = Plots;

# Plot of dose response curves
function gpPlot(ligandName, cellType)
    responseDF = importData()
    time = [0.5, 1, 2, 4]
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[:, 1]

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    #filter!(row -> string(row["Date"]) .== date, filtFrame)

    realDataDF = filtFrame[!, [:Dose, :Time, :Mean]]
    realDataDF = groupby(realDataDF, [:Time, :Dose])
    realDataDF = combine(realDataDF, :Mean => mean)

    fullDataX = filtFrame[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc]]
    intrinsLevels = identity.(convert(Matrix, fullDataX)[1, 3:9])
    xMat = zeros(length(doseVec), length(intrinsLevels) + 2)

    #Train data
    X, y, df = getGPdata()
    gp = gaussianProcess(X', y)
    colors = ["aqua", "coral", "darkorchid", "goldenrod"]
    μs = zeros(length(time), length(doseVec))
    σ²s = similar(μs)

    pl1 = plt.plot()

    for (i, ITtime) in enumerate(time)
        xMat = zeros(length(doseVec), length(intrinsLevels) + 2)
        xMat[:, 1] .= log10.(doseVec)
        xMat[:, 2] .= ITtime
        xMat[:, 3:9] .= repeat(intrinsLevels, outer = [1, length(doseVec)])'
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
        plt.plot!(doseVec, μs[i, :], c = colors[i], xscale = :log10, label = ITtime, legend = :bottomright)

        if length(log10.(realDataDF[realDataDF.Time .== ITtime, :].Mean_mean .+ 1)) > 0
            plt.scatter!(doseVec, log10.(realDataDF[realDataDF.Time .== ITtime, :].Mean_mean .+ 1), c = colors[i], xscale = :log10, label = "")
        end

    end

    ylabel!("pSTAT")
    xlabel!("Dose (nM)")

    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ3()
    l = @layout [a b c d; e f g h; j i k l]
    p1 = gpPlot("WT N-term", "Treg")
    p2 = gpPlot("WT N-term", "Thelper")
    p3 = gpPlot("WT N-term", "NK")
    p4 = gpPlot("WT N-term", "CD8")
    p5 = gpPlot("H16N N-term", "Treg")
    p6 = gpPlot("H16N N-term", "Thelper")
    p7 = gpPlot("H16N N-term", "NK")
    p8 = gpPlot("H16N N-term", "CD8")
    p9 = gpPlot("R38Q N-term", "Treg")
    p10 = gpPlot("R38Q N-term", "Thelper")
    p11 = gpPlot("R38Q N-term", "NK")
    p12 = gpPlot("R38Q N-term", "CD8")
    #draw(SVG("figureJ1.svg", 1000px, 800px), p1)
    ffig = plt.plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, layout = l, size = (1200, 1200))
    plt.savefig(ffig, joinpath(dirname(pathof(gcSolver)), "..", "figureJ3.svg"))
end
