""" This file builds the depletion manuscript, Figure 2. """

const dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")
responseDF = DataFrame!(CSV.File(joinpath(dataDir, "WTMuteinsMoments.csv")))

""" Plot an example isobologram. """
function trialplot2()
    X = [1, 2, 3]
    Y = [1, 2, 3]
    pl = plot(x = X, y = Y)
    return pl
end

# Plot of dose response curves
function doseResPlot2(ligandName, cellType, date, unkVec)
    tps = [0.5, 1, 2, 4] .* 60
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[:, 1]
    receptorDF = getExpression()
    tens = [10, 10, 10, 0, 10]
    cellSpecAbund = receptorDF[!, cellType]
    #realDataDF = DataFrame(Dose = Float64[], time = Float64[], pSTAT = Float64[])
    predictDF = DataFrame(Dose = Float64[], tps = Float64[], sensitivity = Float64[])

    filtFrame = filter(row -> row["Ligand"] .== ligandName, responseDF)
    filter!(row -> row["Cell"] .== cellType, filtFrame)
    filter!(row -> string(row["Date"]) .== date, filtFrame)

    """for ind = 1:size(filtFrame)[1]
        push!(realDataDF, (filtFrame[ind, "Dose"], filtFrame[ind, "Time"], filtFrame[ind, "Mean"]))
    end
    realDataDF = groupby(realDataDF, [:time, :Dose])
    realDataDF = combine(realDataDF, :pSTAT => mean)"""

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
        iterParams = fitParams(doseLevel, unkVec, tens .^ cellSpecAbund, cellType)
        if ligandName != "IL2" && ligandName != "IL15"
            iterParams = mutAffAdjust(iterParams, ligandName)
        end
        #gives you pstat results
        #pstatResults = runCkine(time, iterParams, pSTAT5 = true) .* unkVec[24] .* 1e6
        #println("iterParams = ", iterParams)
        #println("tps = ", tps)
        jacResults = runRecJac(tps, iterParams)
        
        for indx = 1:length(tps)
            #use dataframe and push row into it - enter data into data frame
            push!(predictDF, (dose, tps[indx] / 60, jacResults[indx]))
        end
    end
    
    

    pl1 = plot(
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
    checkInputs(tps, params)
    
    # Sigma is the covariance matrix of the input parameters
    function jacF(x)
        pp = vcat(params[1:24], x, params[30:end])
        return runCkine(tps, pp, pSTAT5 = true)
    end

    jac = zeros(5, length(tps))
    ForwardDiff.jacobian!(jac, jacF, params[25:29])
    println("jac = ", jac)
    return jac
end


    
"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ2()
    fitVec = CSV.read(joinpath(dataDir, "fitTry.csv"))
    fitVec = convert(Vector{Float64}, fitVec[!, :value])
    p1 = doseResPlot2("IL2", "Treg", "2019-03-19", fitVec)
    """p1 = doseResPlot2("IL2", "Treg", "2019-03-19", fitVec)
    p2 = doseResPlot2("IL2", "Thelper", "2019-03-19", fitVec)
    p3 = doseResPlot2("IL2", "NK", "2019-03-15", fitVec)
    p4 = doseResPlot2("IL2", "CD8", "2019-03-15", fitVec)
    p5 = doseResPlot2("WT N-term", "Treg", "2019-04-19", fitVec)
    p6 = doseResPlot2("WT N-term", "Thelper", "2019-04-19", fitVec)
    p7 = doseResPlot2("WT N-term", "NK", "2019-05-02", fitVec)
    p8 = doseResPlot2("WT N-term", "CD8", "2019-05-02", fitVec)
    p9 = doseResPlot2("H16N N-term", "Treg", "2019-04-19", fitVec)
    p10 = doseResPlot2("H16N N-term", "Thelper", "2019-04-19", fitVec)
    p11 = doseResPlot2("H16N N-term", "NK", "2019-05-02", fitVec)
    p12 = doseResPlot2("H16N N-term", "CD8", "2019-05-02", fitVec)
    p13 = doseResPlot2("R38Q N-term", "Treg", "2019-04-19", fitVec)
    p14 = doseResPlot2("R38Q N-term", "Thelper", "2019-04-19", fitVec)
    p15 = doseResPlot2("R38Q N-term", "NK", "2019-05-02", fitVec)
    p16 = doseResPlot2("R38Q N-term", "CD8", "2019-05-02", fitVec)"""
    draw(SVG("figureJ2.svg", 1000px, 800px), p1)
    #draw(SVG("figureJ2.svg", 2000px, 1600px), gridstack([p1 p2 p3 p4; p5 p6 p7 p8; p9 p10 p11 p12; p13 p14 p15 p16]))
end
