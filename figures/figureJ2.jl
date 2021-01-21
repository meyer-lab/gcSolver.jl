""" Cell type specificity plots (ODE model), Figure 2. """

using Gadfly;
using Statistics;
using gcSolver;
using DataFrames;
using ForwardDiff;
using StatsFuns;
import LinearAlgebra: diag;

gdf = Gadfly;

# Plot of dose response curves
function doseResPlot2(ligandName, cellType, date, unkVec, alphaCov = true, biv = false)
    responseDF = gcSolver.importData(false)
    DateFrame = gcSolver.importConvFrame()
    tps = [0.5, 1, 2, 4] .* 60
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[!, :Dose]
    sigma = gcSolver.getSigma(cellType)
    if alphaCov
        predictDF = DataFrame(Dose = Float64[], tps = Float64[], sensitivity = Float64[])
    else
        predictDF = DataFrame(Dose = Float64[], tps = Float64[], Variance = Float64[])
    end

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    filter!(row -> string(row["Date"]) .== date, filtFrame)
    filter!(row -> row["Bivalent"] .== biv, filtFrame)

    if alphaCov
        realDataDF = filtFrame[!, [:Dose, :Time, :Mean, :alphStatCov]]
        realDataDF = groupby(realDataDF, [:Time, :Dose])
        realDataDF = combine(realDataDF, :Mean => mean, :alphStatCov => mean)
    else
        realDataDF = filtFrame[!, [:Dose, :Time, :Mean, :Variance]]
        realDataDF = groupby(realDataDF, [:Time, :Dose])
        realDataDF = combine(realDataDF, :Mean => mean, :Variance => mean)
    end


    for (i, dose) in enumerate(doseVec)
        #check if ligand name is IL2

        if ligandName == "IL15"
            #put into second slot
            doseLevel = [0, dose, 0]
        else
            #put ILdose into first slot
            doseLevel = [dose, 0, 0]
        end

        idxx = findfirst(responseDF.Cell .== cellType)
        iterParams =
            gcSolver.fitParams(doseLevel, unkVec, 10.0 .^ Vector{Float64}(responseDF[idxx, [:IL2Ra, :IL2Rb, :gc, :IL15Ra, :IL7Ra]]), cellType)
        if ligandName != "IL2" && ligandName != "IL15"
            iterParams = gcSolver.mutAffAdjust(iterParams, responseDF[findfirst(responseDF.Ligand .== ligandName), [:IL2RaKD, :IL2RBGKD]])
        end

        if alphaCov
            JacResults = runAlphJac(tps, iterParams, sigma)
            for indx = 1:length(tps)
                #use dataframe and push row into it - enter data into data frame
                push!(predictDF, (dose, tps[indx] / 60, JacResults[indx]))
            end
        else
            VarResults = runCkineVarProp(tps, iterParams, sigma)
            for indx = 1:length(tps)
                #use dataframe and push row into it - enter data into data frame
                push!(predictDF, (dose, tps[indx] / 60, VarResults[indx]))
            end
        end


    end

    if alphaCov
        predictDF.sensitivity .*= filter(row -> row.Date ∈ [date], DateFrame).Conv .^ 2
        pl1 = gdf.plot(
            layer(predictDF, x = :Dose, y = :sensitivity, color = :tps, Geom.line),
            layer(realDataDF, x = :Dose, y = :alphStatCov_mean, color = :Time, Geom.point),
            Scale.x_log10,
            Guide.title(string(cellType, " Response to ", ligandName)),
            Guide.xlabel("Dose"),
            Guide.ylabel("Sensitivity"),
            Scale.color_discrete(),
            Guide.colorkey(title = "Time (hr)", labels = ["4", "2", "1", "0.5"]),
        )
    else
        predictDF.Variance .*= filter(row -> row.Date ∈ [date], DateFrame).Conv .^ 2
        pl1 = gdf.plot(
            layer(predictDF, x = :Dose, y = :Variance, color = :tps, Geom.line),
            layer(realDataDF, x = :Dose, y = :Variance_mean, color = :Time, Geom.point),
            Scale.x_log10,
            Guide.title(string(cellType, " Response to ", ligandName)),
            Guide.xlabel("Dose"),
            Guide.ylabel("Variance"),
            Scale.color_discrete(),
            Guide.colorkey(title = "Time (hr)", labels = ["4", "2", "1", "0.5"]),
            Scale.y_log10,
            Coord.cartesian(ymin = 1, ymax = 10),
        )
    end
    return pl1
end

function runAlphJac(tps::Vector, params::Vector, sigma)
    gcSolver.checkInputs(tps, params)

    # Sigma is the covariance matrix of the input parameters
    function jacF(x)
        pp = vcat(params[1:24], x[1], params[26:end])
        return runCkine(tps, pp, pSTAT5 = true)
    end

    jac = zeros(1, length(tps))
    ForwardDiff.jacobian!(jac, jacF, [params[25]])
    return diag(transpose(jac) * sigma[1, 1] * jac)
end



"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ2()
    fitVec = gcSolver.importFit()
    fitVec = convert(Vector{Float64}, fitVec[!, :Fit])
    fitVec = softplus.(fitVec)

    p1 = doseResPlot2("IL2", "Treg", "3/15/2019", fitVec, true, false)
    p2 = doseResPlot2("IL2", "Thelper", "3/15/2019", fitVec, true, false)
    p3 = doseResPlot2("IL2", "NK", "3/15/2019", fitVec, true, false)
    p4 = doseResPlot2("IL2", "CD8", "3/15/2019", fitVec, true, false)
    p5 = doseResPlot2("N88D C-term", "Treg", "12/5/2019", fitVec, true, false)
    p6 = doseResPlot2("N88D C-term", "Thelper", "12/5/2019", fitVec, true, false)
    p7 = doseResPlot2("N88D C-term", "NK", "12/5/2019", fitVec, true, false)
    p8 = doseResPlot2("N88D C-term", "CD8", "12/5/2019", fitVec, true, false)
    p9 = doseResPlot2("WT N-term", "Treg", "11/8/2019", fitVec, true, false)
    p10 = doseResPlot2("WT N-term", "Thelper", "11/8/2019", fitVec, true, false)
    p11 = doseResPlot2("WT N-term", "NK", "11/8/2019", fitVec, true, false)
    p12 = doseResPlot2("WT N-term", "CD8", "11/8/2019", fitVec, true, false)
    p13 = doseResPlot2("V91K C-term", "Treg", "11/27/2019", fitVec, false, false)
    p14 = doseResPlot2("V91K C-term", "Thelper", "11/27/2019", fitVec, false, false)
    p15 = doseResPlot2("V91K C-term", "NK", "11/27/2019", fitVec, false, false)
    p16 = doseResPlot2("V91K C-term", "CD8", "11/27/2019", fitVec, false, false)
    p17 = doseResPlot2("R38Q N-term", "Treg", "12/5/2019", fitVec, false, false)
    p18 = doseResPlot2("R38Q N-term", "Thelper", "12/5/2019", fitVec, false, false)
    p19 = doseResPlot2("R38Q N-term", "NK", "12/5/2019", fitVec, false, false)
    p20 = doseResPlot2("R38Q N-term", "CD8", "12/5/2019", fitVec, false, false)
    p21 = doseResPlot2("F42Q N-Term", "Treg", "12/5/2019", fitVec, false, false)
    p22 = doseResPlot2("F42Q N-Term", "Thelper", "12/5/2019", fitVec, false, false)
    p23 = doseResPlot2("F42Q N-Term", "NK", "12/5/2019", fitVec, false, false)
    p24 = doseResPlot2("F42Q N-Term", "CD8", "12/5/2019", fitVec, false, false)
    #draw(SVG("figureJ2.svg", 1000px, 800px), p1)
    draw(SVG("figureJ2.svg", 4000px, 2400px), gridstack([p1 p2 p3 p4; p9 p10 p11 p12; p13 p14 p15 p16; p17 p18 p19 p20; p21 p22 p23 p24]))
end
