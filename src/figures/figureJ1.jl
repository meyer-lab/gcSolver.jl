""" This figure builds the experimental/model data overlay. """

using CSV

const dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")
responseDF = CSV.read(joinpath(dataDir, "WTMuteinsMoments.csv"), copycols = true)

""" Plot an example isobologram. """
function trialplot()
    X = [1, 2, 3]
    Y = [1, 2, 3]
    pl = plot(x=X, y=Y);
    return pl
end

# Plot of dose response curves
function doseResPlot(ligandName, cellType, date, unkVec)
    time = [0.5, 1, 2, 4] .* 60
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[:, 1]
    receptorDF = getExpression()
    cellSpecAbund = receptorDF[!, cellType]
    realDataDF = DataFrame(Dose = Float64[], time = Float64[], pSTAT = Float64[])
    predictDF = DataFrame(Dose = Float64[], time = Float64[], pSTAT = Float64[])
    
    for ind = 1:size(responseDF)[1]
        if ligandName == responseDF[ind,"Ligand"]
            if cellType == responseDF[ind,"Cell"]
                if date == string(responseDF[ind,"Date"])   
                    push!(realDataDF, (responseDF[ind,"Dose"], responseDF[ind,"Time"], responseDF[ind,"Mean"]))                    
                end
            end
        end
    end

    for (i, dose) in enumerate(doseVec)
        #check if ligand name is IL2
        if ligandName == "IL2"
            #put ILdose into first slot
            doseLevel = [dose, 0, 0]
        elseif ligandName == "IL15"
            #put into second slot
            doseLevel = [0, dose, 0]
        end
    
        #Gives back 36 parameter long
        iterParams = fitParams(doseLevel, unkVec, cellSpecAbund)
        #gives you pstat results
        pstatResults = runCkine(time, iterParams, pSTAT5 = true) .* unkVec[21] .* 10e6
        for indx = 1:length(time)
        #use dataframe and push row into it - enter data into data frame
            push!(predictDF, (dose, time[indx] / 60, pstatResults[indx]))
        end
    end    
    
    pl1 = plot(layer(realDataDF, x=:Dose, y=:pSTAT, color=:time, Geom.point), layer(predictDF, x=:Dose, y=:pSTAT, color=:time, Geom.line), Scale.x_log10, Scale.y_log10, Guide.title("Dose Response Curves"), Guide.xlabel("Dose"), Guide.ylabel("Pstat Level"), Scale.color_discrete(), Guide.colorkey(title="Time (hr)", labels=["0.5","1","2","4"]))
    
    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ1()
    unkVec = getUnkVec()
    p1 = doseResPlot("IL2", "Thelper", "2019-03-19", unkVec)
    draw(SVG("figureJ1.svg", 1000px, 800px), p1)
    
    #p1 = trialplot()
    #p2 = trialplot()
    #p3 = trialplot()
    #p4 = trialplot()
    

    #draw(SVG("figureJ1.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end