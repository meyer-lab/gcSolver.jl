""" This file builds the depletion manuscript, Figure 4. """

function cellTypeContr(gp, realType, compType, recExp=false)
    x, y, df = getGPdata()
    predDF = DataFrame(realPred = Float64[], fakePred = Float64[], ligand = String[]) 
    
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
    
    ligands = df[df.Cell .== realType, :].Ligand
    #println("ligands = ", ligands)

  
    """#assuming first tuple returned by predict f is the inputs and the second tuple returned is the predictions
    realVals = realPreds[2]
    fakeVals = compPreds[2]"""
    
    # first tuple returned by predict_f is the predictions and the second tuple returned is the standard deviations
    realVals = realPreds[1]
    fakeVals = compPreds[1]
    #for i in enumerate(realVals)
    for i = 1:length(ligands)
        #push!(predDF, (realPreds[i], compPreds[i]), ::float64)
        real = realVals[i]
        fake = fakeVals[i]
        lig = ligands[i]
        push!(predDF, (real, fake, lig))

    end
    #println("predDF = ", predDF[1:20,:])
    predComp = gdf.plot(
        layer(predDF, x = :realPred, y = :fakePred, color = :ligand, Geom.point),
        Guide.title(string(compType, " influence on ", realType, " prediction")),
        Guide.xlabel("Correct Pred"),
        Guide.ylabel("Incorrect Pred"),
        Scale.color_discrete(),
        Guide.colorkey(title = "Muteins", labels = ["IL2", "IL15", "WT N-term", "R38Q N-term", "H16N N-term", "R38Q/H16N"])
    )
        
    
    return predComp
end

"""Generate plots"""
function figureJ4()
    x, y, df = getGPdata()
    trainedGP = gaussianProcess(x', y)
    
    p1 = cellTypeContr(trainedGP, "Treg", "Thelper", false)
    p2 = cellTypeContr(trainedGP, "Treg", "NK", false)
    p3 = cellTypeContr(trainedGP, "Treg", "CD8", false)
    
    p4 = cellTypeContr(trainedGP, "Thelper", "Treg", false)
    p5 = cellTypeContr(trainedGP, "Thelper", "NK", false)
    p6 = cellTypeContr(trainedGP, "Thelper", "CD8", false)
    
    p7 = cellTypeContr(trainedGP, "NK", "Treg", false)
    p8 = cellTypeContr(trainedGP, "NK", "Thelper", false)
    p9 = cellTypeContr(trainedGP, "NK", "CD8", false)
    
    p10 = cellTypeContr(trainedGP, "CD8", "Treg", false)
    p11 = cellTypeContr(trainedGP, "CD8", "Thelper", false)
    p12 = cellTypeContr(trainedGP, "CD8", "NK", false)
    
    #draw(SVG("figureJ4.svg", 4000px, 1600px), p1)
    #draw(SVG("figureJ2.svg", 4000px, 1600px), gridstack([p13 p14 p15 p16; p17 p18 p19 p20; p21 p22 p23 p24]))
    
    draw(SVG("figureJ4.svg", 2500px, 1600px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p11 p12]))
end
