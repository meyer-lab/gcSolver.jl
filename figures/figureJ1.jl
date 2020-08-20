""" This file builds the depletion manuscript, Figure 1. """

using DataFrames;
using gcSolver;
using Gadfly;
using Statistics
gdf = Gadfly;

# Plot of dose response curves
function doseResPlot(ligandName, cellType, date, unkVec)
    responseDF = gcSolver.importData()
    time = [0.5, 1, 2, 4] .* 60
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

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
        iterParams =
            gcSolver.fitParams(doseLevel, unkVec, 10.0 .^ Vector{Float64}(responseDF[idxx, [:IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc]]), cellType)
        if ligandName != "IL2" && ligandName != "IL15"
            iterParams = gcSolver.mutAffAdjust(iterParams, responseDF[findfirst(responseDF.Ligand .== ligandName), [:IL2RaKD, :IL2RBGKD]])
        end
        #gives you pstat results
        pstatResults = gcSolver.runCkine(time, iterParams, pSTAT5 = true) .* unkVec[24] .* 1e6
        for indx = 1:length(time)
            #use dataframe and push row into it - enter data into data frame
            push!(predictDF, (dose, time[indx] / 60, pstatResults[indx]))
        end
    end

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

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ1()
    fitVec = gcSolver.importFit()
    fitVec = convert(Vector{Float64}, fitVec[!, :Fit])
    """
    p1 = doseResPlot("IL2", "Treg", "3/19/2019", fitVec)
    p2 = doseResPlot("IL2", "Thelper", "3/19/2019", fitVec)
    p3 = doseResPlot("IL2", "NK", "3/15/2019", fitVec)
    p4 = doseResPlot("IL2", "CD8", "3/15/2019", fitVec)
    p5 = doseResPlot("IL15", "Treg", "3/19/2019", fitVec)
    p6 = doseResPlot("IL15", "Thelper", "3/19/2019", fitVec)
    p7 = doseResPlot("IL15", "NK", "3/15/2019", fitVec)
    p8 = doseResPlot("IL15", "CD8", "3/15/2019", fitVec)
    """
    p9 = doseResPlot("R38Q/H16N", "Treg", "4/19/2019", fitVec)
    p10 = doseResPlot("R38Q/H16N", "Thelper", "4/19/2019", fitVec)
    p11 = doseResPlot("R38Q/H16N", "NK", "5/02/2019", fitVec)
    p12 = doseResPlot("R38Q/H16N", "CD8", "3/19/2019", fitVec)
    
    p13 = doseResPlot("WT N-term", "Treg", "4/19/2019", fitVec)
    p14 = doseResPlot("WT N-term", "Thelper", "4/19/2019", fitVec)
    p15 = doseResPlot("WT N-term", "NK", "5/2/2019", fitVec)
    p16 = doseResPlot("WT N-term", "CD8", "5/2/2019", fitVec)
    p17 = doseResPlot("H16N N-term", "Treg", "4/19/2019", fitVec)
    p18 = doseResPlot("H16N N-term", "Thelper", "4/19/2019", fitVec)
    p19 = doseResPlot("H16N N-term", "NK", "5/2/2019", fitVec)
    p20 = doseResPlot("H16N N-term", "CD8", "5/2/2019", fitVec)
    p21 = doseResPlot("R38Q N-term", "Treg", "4/19/2019", fitVec)
    p22 = doseResPlot("R38Q N-term", "Thelper", "4/19/2019", fitVec)
    p23 = doseResPlot("R38Q N-term", "NK", "5/2/2019", fitVec)
    p24 = doseResPlot("R38Q N-term", "CD8", "5/2/2019", fitVec)
    #draw(SVG("figureJ2.svg", 1000px, 800px), p1)
    draw(SVG("figureJ1.svg", 4000px, 1600px), gridstack([p9 p10 p11 p12; p13 p14 p15 p16; p17 p18 p19 p20; p21 p22 p23 p24]))
end
