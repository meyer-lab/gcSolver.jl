""" This file builds the depletion manuscript, Figure 2. """

using Gadfly;
using Statistics;
using gcSolver;
using DataFrames;
using ForwardDiff;
gdf = Gadfly;

# Plot of dose response curves
function doseResPlot2(ligandName, cellType, date, unkVec)
    responseDF = gcSolver.importData(true)
    tps = [0.5, 1, 2, 4] .* 60
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]

    predictDF = DataFrame(Dose = Float64[], tps = Float64[], sensitivity = Float64[])

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    filter!(row -> string(row["Date"]) .== date, filtFrame)

    realDataDF = filtFrame[!, [:Dose, :Time, :Mean]]
    realDataDF = groupby(realDataDF, [:Time, :Dose])
    realDataDF = combine(realDataDF, :Mean => mean)


    for (i, dose) in enumerate(doseVec)
        #check if ligand name is IL2

        if ligandName == "IL15"
            #put into second slot
            doseLevel = [0, dose, 0]
        else
            #put ILdose into first slot
            doseLevel = [dose, 0, 0]
        end

        #Gives back 36 parameter long
        idxx = findfirst(responseDF.Cell .== cellType)
        iterParams =
            gcSolver.fitParams(doseLevel, unkVec, 10.0 .^ Vector{Float64}(responseDF[idxx, [:IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc]]), cellType)
        if ligandName != "IL2" && ligandName != "IL15"
            iterParams = gcSolver.mutAffAdjust(iterParams, responseDF[findfirst(responseDF.Ligand .== ligandName), [:IL2RaKD, :IL2RBGKD]])
        end

        jacResults = runRecJac(tps, iterParams)

        for indx = 1:length(tps)
            #use dataframe and push row into it - enter data into data frame
            push!(predictDF, (dose, tps[indx] / 60, jacResults[1, indx]))
        end
    end

    pl1 = gdf.plot(
        layer(predictDF, x = :Dose, y = :sensitivity, color = :tps, Geom.line),
        Scale.x_log10,
        Guide.title(string(cellType, " Response to ", ligandName)),
        Guide.xlabel("Dose"),
        Guide.ylabel("Sensitivity"),
        Scale.color_discrete(),
        Guide.colorkey(title = "Time (hr)", labels = ["4", "2", "1", "0.5"]),
    )

    return pl1
end

function runRecJac(tps::Vector, params::Vector)
    gcSolver.checkInputs(tps, params)

    # Sigma is the covariance matrix of the input parameters
    function jacF(x)
        pp = vcat(params[1:24], x, params[30:end])
        return runCkine(tps, pp, pSTAT5 = true)
    end

    jac = zeros(5, length(tps))
    ForwardDiff.jacobian!(jac, jacF, params[25:29])
    return jac
end



"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ2()
    fitVec = gcSolver.importFit()
    fitVec = convert(Vector{Float64}, fitVec[!, :Fit])
    fitVec = softplus.(fitVec)
    
    p1 = doseResPlot2("IL2", "Treg", "3/19/2019", fitVec)
    p2 = doseResPlot2("IL2", "Thelper", "3/19/2019", fitVec)
    p3 = doseResPlot2("IL2", "NK", "3/15/2019", fitVec)
    p4 = doseResPlot2("IL2", "CD8", "3/15/2019", fitVec)
    p5 = doseResPlot2("IL15", "Treg", "3/19/2019", fitVec)
    p6 = doseResPlot2("IL15", "Thelper", "3/19/2019", fitVec)
    p7 = doseResPlot2("IL15", "NK", "3/15/2019", fitVec)
    p8 = doseResPlot2("IL15", "CD8", "3/15/2019", fitVec)
    p9 = doseResPlot2("WT C-term", "Treg", "3/1/2019", fitVec)
    p10 = doseResPlot2("WT C-term", "Thelper", "3/1/2019", fitVec)
    p11 = doseResPlot2("WT C-term", "NK", "3/1/2019", fitVec)
    p12 = doseResPlot2("WT C-term", "CD8", "3/1/2019", fitVec)
    p13 = doseResPlot2("WT N-term", "Treg", "3/1/2019", fitVec)
    p14 = doseResPlot2("WT N-term", "Thelper", "3/1/2019", fitVec)
    p15 = doseResPlot2("WT N-term", "NK", "3/1/2019", fitVec)
    p16 = doseResPlot2("WT N-term", "CD8", "3/1/2019", fitVec)
    p17 = doseResPlot2("V91K C-term", "Treg", "3/1/2019", fitVec)
    p18 = doseResPlot2("V91K C-term", "Thelper", "3/1/2019", fitVec)
    p19 = doseResPlot2("V91K C-term", "NK", "3/1/2019", fitVec)
    p20 = doseResPlot2("V91K C-term", "CD8", "3/1/2019", fitVec)
    p21 = doseResPlot2("R38Q N-term", "Treg", "3/1/2019", fitVec)
    p22 = doseResPlot2("R38Q N-term", "Thelper", "3/1/2019", fitVec)
    p23 = doseResPlot2("R38Q N-term", "NK", "3/1/2019", fitVec)
    p24 = doseResPlot2("R38Q N-term", "CD8", "3/1/2019", fitVec)
    p25 = doseResPlot2("F42Q N-Term", "Treg", "3/1/2019", fitVec)
    p26 = doseResPlot2("F42Q N-Term", "Thelper", "3/1/2019", fitVec)
    p27 = doseResPlot2("F42Q N-Term", "NK", "3/1/2019", fitVec)
    p28 = doseResPlot2("F42Q N-Term", "CD8", "3/1/2019", fitVec)
    p29 = doseResPlot2("N88D C-term", "Treg", "3/1/2019", fitVec)
    p30 = doseResPlot2("N88D C-term", "Thelper", "3/1/2019", fitVec)
    p31 = doseResPlot2("N88D C-term", "NK", "3/1/2019", fitVec)
    p32 = doseResPlot2("N88D C-term", "CD8", "3/1/2019", fitVec)
    #draw(SVG("figureJ2.svg", 1000px, 800px), p1)
    draw(SVG("figureJ2.svg", 4000px, 2400px), gridstack([p1 p2 p3 p4; p5 p6 p7 p8; p9 p10 p11 p12; p13 p14 p15 p16; p17 p18 p19 p20; p21 p22 p23 p24; p25 p26 p27 p28; p29 p30 p31 p32]))
end
