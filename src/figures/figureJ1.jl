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
function doseResPlot(ligandName, cellType, date)
    Xhalf = zeros(12)
    Yhalf = zeros(12)
    X1 = zeros(12)
    Y1 = zeros(12)
    X2 = zeros(12)
    Y2 = zeros(12)
    X4 = zeros(12)
    Y4 = zeros(12)
    iHalf = 1
    i1 = 1
    i2 = 1
    i4 = 1
    for ind = 1:2640
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
    
    pl1 = plot(layer(x=Xhalf, y=Yhalf, Theme(default_color="green")), layer(x=X1, y=Y1, Theme(default_color="blue")), layer(x=X2, y=Y2, Theme(default_color="red")), layer(x=X4, y=Y4, Theme(default_color="orange")), Guide.manual_color_key("Legend", ["0.5 Hours", "1 Hour", "2 Hours", "4 Hours"], ["green", "blue", "red", "orange"]), Guide.title("Dose Response Curves at Various Times"), Guide.xlabel("Dose"), Guide.ylabel("Pstat Level"))
    
    return pl1
end

"""Use this if you want to change the parameters here and not input any in the command line"""
function figureJ1()
    p1 = doseResPlot("IL2", "Thelper", "2019-03-19")
    draw(SVG("figureJ1.svg", 1000px, 800px), p1)
    
    """p1 = trialplot()
    p2 = trialplot()
    p3 = trialplot()
    p4 = trialplot()
    

    draw(SVG("figureJ1.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))"""
end

"""Use this if you want to input parameters in command line"""
function doseResFig(ligandName, cellType, date)
    p1 = doseResPlot(ligandName, cellType, date)
    draw(SVG("figureJ1.svg", 1000px, 800px), p1)
end