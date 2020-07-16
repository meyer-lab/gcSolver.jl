""" This file builds the depletion manuscript, Figure 4. """

function cellTypeContr(gp, realType, compType, recExp=false)
    x, y, df = getGPdata()
    predDF = DataFrame(realPred = Float64[], fakePred = Float64[])
    
    
    #RealX = only rows with x axis cell-type, want it to make predictions
    realX = x[df.Cell .== realType, :]    
    #realX = filter(row -> row["Cell"] .!== realType, df)
    
    #println("RealX size = ", size(realX,2))
    
    realPreds = predict_f(gp, realX')
    # returns tuple of arrays
    
    #modify RealX
    hotEnc = cellHotEnc(compType)
    realX[:, 10:13] = repeat(hotEnc, outer = [1, size(realX)[1]])'

    #realX[10:13,:] = cellHotEnc(compType)
    compPreds = predict_f(gp, realX')

  
    """#assuming first tuple returned by predict f is the inputs and the second tuple returned is the predictions
    realVals = realPreds[2]
    fakeVals = compPreds[2]"""
    
    # first tuple returned by predict_f is the predictions and the second tuple returned is the standard deviations
    realVals = realPreds[1]
    fakeVals = compPreds[1]
    #for i in enumerate(realVals)
    for i = 1:length(realVals)
        #push!(predDF, (realPreds[i], compPreds[i]), ::float64)
        real = realVals[i]
        fake = fakeVals[i]
        push!(predDF, (real, fake))
    end
    
    #println("predDF = ", predDF[1:20,:])
    predComp = gdf.plot(
        layer(predDF, x = :realPred, y = :fakePred, Geom.point),
        Guide.title(string(compType, " influence on ", realType, " prediction")),
        Guide.xlabel("Correct Pred"),
        Guide.ylabel("Incorrect Pred")
    )
        
    
    return predComp
end

"""Generate plots"""
function figureJ4()
    x, y, df = getGPdata()
    trainedGP = gaussianProcess(x', y)
    
    p1 = cellTypeContr(trainedGP, "Treg", "Thelper", true)
    p2 = cellTypeContr(trainedGP, "Treg", "NK")
    p3 = cellTypeContr(trainedGP, "Treg", "CD8")
    
    p4 = cellTypeContr(trainedGP, "Thelper", "Treg")
    p5 = cellTypeContr(trainedGP, "Thelper", "NK")
    p6 = cellTypeContr(trainedGP, "Thelper", "CD8")
    
    p7 = cellTypeContr(trainedGP, "NK", "Treg")
    p8 = cellTypeContr(trainedGP, "NK", "Thelper")
    p9 = cellTypeContr(trainedGP, "NK", "CD8")
    
    p10 = cellTypeContr(trainedGP, "CD8", "Treg")
    p11 = cellTypeContr(trainedGP, "CD8", "Thelper")
    p12 = cellTypeContr(trainedGP, "CD8", "NK")
    
    draw(SVG("figureJ4.svg", 2500px, 1600px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p11 p12]))
end
