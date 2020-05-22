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

""" Plot of dose response curves """
function doseResPlot(ligandName, cellType, date, unkVec)
    time = [0.5, 1, 2, 4] .* 60
    doseVec = unique(responseDF, "Dose")
    doseVec = doseVec[:, 1]
    receptorDF = getExpression()
    """cellSpecAbund = receptorDF[!, symbol(cellType)] that Brian gave me says symbol not defined""" 
    cellSpecAbund = receptorDF[!, cellType]
    """predictDF = DataFrames(Dose = Float64[], time = Float64[], pSTAT = Float64[]) says objects of type Module are not callable"""
    predictDF = DataFrame(Dose = Float64[], time = Float64[], pSTAT = Float64[])
    Xhalf = zeros(12)
    Yhalf = copy(Xhalf)
    X1 = copy(Xhalf)
    Y1 = copy(Xhalf)
    X2 = copy(Xhalf)
    Y2 = copy(Xhalf)
    X4 = copy(Xhalf)
    Y4 = copy(Xhalf)
    iHalf = 1
    i1 = 1
    i2 = 1
    i4 = 1
    for ind = 1:size(responseDF)[1]
        ligand = responseDF[ind,"Ligand"]
        if ligand == ligandName
            cell = responseDF[ind,"Cell"]
            if cell == cellType
                """println("cellsuccess")"""
                day = responseDF[ind,"Date"]
                """println(typeof(day))"""
                """println(typeof(date))"""
                dayStr = string(day)
                """println(typeof(dayStr))"""
                if dayStr == date
                    """println("datesuccess")"""
                    if responseDF[ind,"Time"] == 0.5
                        """println("timesuccess")"""
                        Xhalf[iHalf] = responseDF[ind,"Dose"]
                        Yhalf[iHalf] = responseDF[ind,"Mean"]
                        iHalf += 1;
                    end
                    if responseDF[ind,"Time"] == 1
                        X1[i1] = responseDF[ind,"Dose"]
                        Y1[i1] = responseDF[ind,"Mean"]
                        i1 += 1;
                    end
                    if responseDF[ind,"Time"] == 2
                        X2[i2] = responseDF[ind,"Dose"]
                        Y2[i2] = responseDF[ind,"Mean"]
                        i2 += 1;
                    end
                    if responseDF[ind,"Time"] == 4
                        X4[i4] = responseDF[ind,"Dose"]
                        Y4[i4] = responseDF[ind,"Mean"]
                        i4 += 1;
                    end
                end
            end
        end
    end
    
    #for ILdose in doseVec gave error - AbstractDataFrame is not iterable. Use eachrow(df) to get a row iterator or eachcol(df) to get a column iterator

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
        pstatResults = runCkine(time, iterParams, pSTAT5 = true)
        for indx = 1:length(time)
        #use dataframe and push row into it - enter data into data frame
            push!(predictDF, (dose, time[indx], pstatResults[indx]))
        end
    end
    #should contain all 48 of prediction^, next step is plot it lines using gadfly (4 separate lines depending on time entry), should be something like plot(layer = line, data = predictDF, x=dose, y = pstatresponse, hue/color=Time), prob using predictDF
    
            
    pl1 = plot(layer(x=Xhalf, y=Yhalf, Theme(default_color="green")), layer(x=X1, y=Y1, Theme(default_color="blue")), layer(x=X2, y=Y2, Theme(default_color="red")), layer(x=X4, y=Y4, Theme(default_color="orange")), Guide.manual_color_key("Legend", ["0.5 Hours", "1 Hour", "2 Hours", "4 Hours"], ["green", "blue", "red", "orange"]), Guide.title("Dose Response Curves at Various Times"), Guide.xlabel("Dose"), Guide.ylabel("Pstat Level"))
    
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