""" This file builds the depletion manuscript, Figure 1. """

# Plot of dose response curves
function doseResPlot(ligandName, cellType, date, unkVec)
    responseDF = importData()
    time = [0.5, 1, 2, 4] .* 60
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[:, 1]

    predictDF = DataFrame(Dose = Float64[], time = Float64[], pSTAT = Float64[])

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    filter!(row -> string(row["Date"]) .== date, filtFrame)

    realDataDF = filtFrame[!, [:Dose, :Time, :Mean]]
    realDataDF = groupby(realDataDF, [:Time, :Dose])
    realDataDF = combine(realDataDF, :Mean => mean)

    for (i, dose) in enumerate(doseVec)

        # Place ligand dose
        if ligandName == "IL15"
            doseLevel = [0, dose, 0]
        else
            doseLevel = [dose, 0, 0]
        end

        #Gives back 36 parameter long
        idxx = findfirst(responseDF.Cell .== cellType)
        iterParams = fitParams(doseLevel, unkVec, 10.0 .^ Vector{Float64}(responseDF[idxx, [:IL15Ra, :IL2Ra , :IL2Rb , :IL7Ra, :gc]]), cellType)
        if ligandName != "IL2" && ligandName != "IL15"
            iterParams = mutAffAdjust(iterParams, ligandName)
        end
        #gives you pstat results
        pstatResults = runCkine(time, iterParams, pSTAT5 = true) .* unkVec[24] .* 1e6
        for indx = 1:length(time)
            #use dataframe and push row into it - enter data into data frame
            push!(predictDF, (dose, time[indx] / 60, pstatResults[indx]))
        end
    end

    pl1 = plot(
        layer(realDataDF, x = :Dose, y = :Mean_mean, color = :Time, Geom.point),
        layer(predictDF, x = :Dose, y = :pSTAT, color = :time, Geom.line),
        Scale.x_log10,
        Guide.title(string(cellType, " Response to ", ligandName)),
        Guide.xlabel("Dose"),
        Guide.ylabel("pSTAT Level"),
        Scale.color_discrete(),
        Guide.colorkey(title = "Time (hr)", labels = ["4", "2", "1", "0.5"]),
    )

    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ1()
    fitVec = importFit()
    fitVec = convert(Vector{Float64}, fitVec[!, :value])
    p1 = doseResPlot("IL2", "Treg", "2019-03-19", fitVec)
    p2 = doseResPlot("IL2", "Thelper", "2019-03-19", fitVec)
    p3 = doseResPlot("IL2", "NK", "2019-03-15", fitVec)
    p4 = doseResPlot("IL2", "CD8", "2019-03-15", fitVec)
    p5 = doseResPlot("WT N-term", "Treg", "2019-04-19", fitVec)
    p6 = doseResPlot("WT N-term", "Thelper", "2019-04-19", fitVec)
    p7 = doseResPlot("WT N-term", "NK", "2019-05-02", fitVec)
    p8 = doseResPlot("WT N-term", "CD8", "2019-05-02", fitVec)
    p9 = doseResPlot("H16N N-term", "Treg", "2019-04-19", fitVec)
    p10 = doseResPlot("H16N N-term", "Thelper", "2019-04-19", fitVec)
    p11 = doseResPlot("H16N N-term", "NK", "2019-05-02", fitVec)
    p12 = doseResPlot("H16N N-term", "CD8", "2019-05-02", fitVec)
    p13 = doseResPlot("R38Q N-term", "Treg", "2019-04-19", fitVec)
    p14 = doseResPlot("R38Q N-term", "Thelper", "2019-04-19", fitVec)
    p15 = doseResPlot("R38Q N-term", "NK", "2019-05-02", fitVec)
    p16 = doseResPlot("R38Q N-term", "CD8", "2019-05-02", fitVec)
    #draw(SVG("figureJ1.svg", 1000px, 800px), p1)
    draw(SVG("figureJ1.svg", 2000px, 1600px), gridstack([p1 p2 p3 p4; p5 p6 p7 p8; p9 p10 p11 p12; p13 p14 p15 p16]))
end
