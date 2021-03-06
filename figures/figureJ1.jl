""" Dose response curve (ODE model), Figure 1. """

using DataFrames;
using gcSolver;
using Gadfly;
using Statistics;
using StatsFuns;
gdf = Gadfly;

# Plot of dose response curves
function doseResPlot(ligandName, cellType, date, unkVec, biv = false)
    responseDF = gcSolver.importData(false)
    time = [0.5, 1, 2, 4] .* 60
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

    predictDF = DataFrame(Dose = Float64[], time = Float64[], pSTAT = Float64[])

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    filter!(row -> row["Date"] .== date, filtFrame)
    filter!(row -> row["Bivalent"] .== biv, filtFrame)


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
        iterParams =
            gcSolver.fitParams(doseLevel, unkVec, 10.0 .^ Vector{Float64}(responseDF[idxx, [:IL2Ra, :IL2Rb, :gc, :IL15Ra, :IL7Ra]]), cellType)
        if ligandName != "IL2" && ligandName != "IL15"
            iterParams = gcSolver.mutAffAdjust(iterParams, responseDF[findfirst(responseDF.Ligand .== ligandName), [:IL2RaKD, :IL2RBGKD]])
        end
        #gives you pstat results
        pstatResults = gcSolver.runCkine(time, iterParams, pSTAT5 = true)
        for indx = 1:length(time)
            #use dataframe and push row into it - enter data into data frame
            push!(predictDF, (dose, time[indx] / 60, pstatResults[indx]))
        end
    end

    DateFrame = gcSolver.importConvFrame()
    predictDF.pSTAT .*= filter(row -> row.Date ∈ [date], DateFrame).Conv

    pl1 = gdf.plot(
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


"""This function compares the fit of Farhat and our Jak/STAT pathway model"""
function fitCompare(fitVec)
    numParams = 18
    farhatParams = gcSolver.getFarhatVec()[1:numParams]
    fitVecComp = fitVec[1:numParams]
    paramNames = gcSolver.getParamNames()[1:numParams]
    compDF = DataFrame(
        Parameter = repeat(paramNames, 2),
        Value = append!(farhatParams, fitVecComp),
        Source = append!(fill("Farhat", numParams), fill("Jak/STAT", numParams)),
    )
    pl1 = gdf.plot(compDF, y = :Value, x = :Parameter, color = :Source, Scale.y_log10)
    return pl1
end


"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ1()
    fitVec = gcSolver.importFit()
    fitVec = convert(Vector{Float64}, fitVec[!, :Fit])
    fitVec = softplus.(fitVec)

    pParam = fitCompare(fitVec)
    p1 = doseResPlot("IL2", "Treg", "4/18/2019", fitVec, 0)
    p2 = doseResPlot("IL2", "Thelper", "4/18/2019", fitVec, 0)
    p3 = doseResPlot("IL2", "NK", "4/18/2019", fitVec, 0)
    p4 = doseResPlot("IL2", "CD8", "4/18/2019", fitVec, 0)
    p5 = doseResPlot("N88D C-term", "Treg", "12/5/2019", fitVec, 0)
    p6 = doseResPlot("N88D C-term", "Thelper", "12/5/2019", fitVec, 0)
    p7 = doseResPlot("N88D C-term", "NK", "12/5/2019", fitVec, 0)
    p8 = doseResPlot("N88D C-term", "CD8", "12/5/2019", fitVec, 0)
    p9 = doseResPlot("WT N-term", "Treg", "11/8/2019", fitVec, 0)
    p10 = doseResPlot("WT N-term", "Thelper", "11/8/2019", fitVec, 0)
    p11 = doseResPlot("WT N-term", "NK", "11/8/2019", fitVec, 0)
    p12 = doseResPlot("WT N-term", "CD8", "11/8/2019", fitVec, 0)
    p13 = doseResPlot("V91K C-term", "Treg", "11/27/2019", fitVec, 0)
    p14 = doseResPlot("V91K C-term", "Thelper", "11/27/2019", fitVec, 0)
    p15 = doseResPlot("V91K C-term", "NK", "11/27/2019", fitVec, 0)
    p16 = doseResPlot("V91K C-term", "CD8", "11/27/2019", fitVec, 0)
    p17 = doseResPlot("R38Q N-term", "Treg", "12/5/2019", fitVec, 0)
    p18 = doseResPlot("R38Q N-term", "Thelper", "12/5/2019", fitVec, 0)
    p19 = doseResPlot("R38Q N-term", "NK", "12/5/2019", fitVec, 0)
    p20 = doseResPlot("R38Q N-term", "CD8", "12/5/2019", fitVec, 0)
    p21 = doseResPlot("F42Q N-Term", "Treg", "12/5/2019", fitVec, 0)
    p22 = doseResPlot("F42Q N-Term", "Thelper", "12/5/2019", fitVec, 0)
    p23 = doseResPlot("F42Q N-Term", "NK", "12/5/2019", fitVec, 0)
    p24 = doseResPlot("F42Q N-Term", "CD8", "12/5/2019", fitVec, 0)
    p25 = doseResPlot("IL15", "Treg", "4/18/2019", fitVec, 0)
    p26 = doseResPlot("IL15", "Thelper", "4/18/2019", fitVec, 0)
    p27 = doseResPlot("IL15", "NK", "4/18/2019", fitVec, 0)
    p28 = doseResPlot("IL15", "CD8", "4/18/2019", fitVec, 0)

    #draw(SVG("figureJ2.svg", 1000px, 800px), p1)
    draw(
        SVG("figureJ1.svg", 4000px, 2400px),
        gridstack(
            [
                pParam pParam pParam pParam
                p1 p2 p3 p4
                p5 p6 p7 p8
                p9 p10 p11 p12
                p13 p14 p15 p16
                p17 p18 p19 p20
                p21 p22 p23 p24
                p25 p26 p27 p28
            ],
        ),
    )
end
