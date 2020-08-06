""" This file builds the depletion manuscript, Figure 4. """

using Gadfly;
using gcSolver;
using DataFrames;
using GaussianProcesses;
gdf = Gadfly;

function cellTypeContr(gp, realType, compType, recExp = false)

    x, y, df = gcSolver.getGPdata()

    predDF = DataFrame(realPred = Float64[], fakePred = Float64[], ligand = String[])

    #RealX = only rows with x axis cell-type, want it to make predictions
    realX = x[df.Cell .== realType, :]
    #realX = filter(row -> row["Cell"] .!== realType, df)

    #println("RealX size = ", size(realX,2))

    realPreds = predict_f(gp, realX')
    # returns tuple of arrays

    if recExp == true
        compExpression = x[df.Cell .== compType, :][1, 5:9]
        realX[:, 5:9] = repeat(compExpression, outer = [1, size(realX)[1]])'
        """elseif CD4 == false
           #modify RealX
            hotEnc = cellHotEnc(compType)
            realX[:, 10:13] = repeat(hotEnc, outer = [1, size(realX)[1]])'
            #realX[10:13,:] = cellHotEnc(compType)
        elseif CD4 == true"""
    else
        hotEnc = gcSolver.cellHotEnc(compType)
        """if compType == "Treg"
            hotEnc = cellHotEnc("Thelper")
        end"""
        realX[:, 10:size(realX, 2)] = repeat(hotEnc, outer = [1, size(realX)[1]])'
        #println("size (cols) = ", size(realX,2))
    end

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

    if recExp == true
        predComp = gdf.plot(
            layer(predDF, x = :realPred, y = :fakePred, color = :ligand, Geom.point),
            Guide.title(string(compType, " influence on ", realType, " prediction (Receptor profile, CD4+)")),
            Guide.xlabel("Correct Pred"),
            Guide.ylabel("Incorrect Pred"),
            Scale.color_discrete(),
            Guide.colorkey(title = "Muteins"),
        )
    else
        predComp = gdf.plot(
            layer(predDF, x = :realPred, y = :fakePred, color = :ligand, Geom.point),
            Guide.title(string(compType, " influence on ", realType, " prediction (CD4+)")),
            Guide.xlabel("Correct Pred"),
            Guide.ylabel("Incorrect Pred"),
            Scale.color_discrete(),
            Guide.colorkey(title = "Muteins"),
        )
    end



    return predComp
end

"""Generate plots"""
function figureJ4()
    x, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(x', y)

    #x4, y4, df4 = getGPdata(true)
    #trainedGP4 = gaussianProcess(x4', y4)

    p1 = cellTypeContr(trainedGP, "Treg", "Thelper", false)
    p2 = cellTypeContr(trainedGP, "Treg", "NK", false)
    p3 = cellTypeContr(trainedGP, "Treg", "CD8", false)
    p1rec = cellTypeContr(trainedGP, "Treg", "Thelper", true)
    p2rec = cellTypeContr(trainedGP, "Treg", "NK", true)
    p3rec = cellTypeContr(trainedGP, "Treg", "CD8", true)
    """p1CD4 = cellTypeContr(trainedGP4, "Treg", "Thelper", false, true)
    p2CD4 = cellTypeContr(trainedGP4, "Treg", "NK", false, true)
    p3CD4 = cellTypeContr(trainedGP4, "Treg", "CD8", false, true)
    p1CD4rec = cellTypeContr(trainedGP4, "Treg", "Thelper", true, true)
    p2CD4rec = cellTypeContr(trainedGP4, "Treg", "NK", true, true)
    p3CD4rec = cellTypeContr(trainedGP4, "Treg", "CD8", true, true)"""

    p4 = cellTypeContr(trainedGP, "Thelper", "Treg")
    p5 = cellTypeContr(trainedGP, "Thelper", "NK")
    p6 = cellTypeContr(trainedGP, "Thelper", "CD8")
    p4rec = cellTypeContr(trainedGP, "Thelper", "Treg", true)
    p5rec = cellTypeContr(trainedGP, "Thelper", "NK", true)
    p6rec = cellTypeContr(trainedGP, "Thelper", "CD8", true)
    """p4CD4 = cellTypeContr(trainedGP4, "Thelper", "Treg", false, true)
    p5CD4 = cellTypeContr(trainedGP4, "Thelper", "NK", false, true)
    p6CD4 = cellTypeContr(trainedGP4, "Thelper", "CD8", false, true)
    p4CD4rec = cellTypeContr(trainedGP4, "Thelper", "Treg", true, true)
    p5CD4rec = cellTypeContr(trainedGP4, "Thelper", "NK", true, true)
    p6CD4rec = cellTypeContr(trainedGP4, "Thelper", "CD8", true, true)"""

    p7 = cellTypeContr(trainedGP, "NK", "Treg")
    p8 = cellTypeContr(trainedGP, "NK", "Thelper")
    p9 = cellTypeContr(trainedGP, "NK", "CD8")
    p7rec = cellTypeContr(trainedGP, "NK", "Treg", true)
    p8rec = cellTypeContr(trainedGP, "NK", "Thelper", true)
    p9rec = cellTypeContr(trainedGP, "NK", "CD8", true)
    """p7CD4 = cellTypeContr(trainedGP4, "NK", "Treg", false, true)
    p8CD4 = cellTypeContr(trainedGP4, "NK", "Thelper", false, true)
    p9CD4 = cellTypeContr(trainedGP4, "NK", "CD8", false, true)
    p7CD4rec = cellTypeContr(trainedGP4, "NK", "Treg", true, true)
    p8CD4rec = cellTypeContr(trainedGP4, "NK", "Thelper", true, true)
    p9CD4rec = cellTypeContr(trainedGP4, "NK", "CD8", true, true)"""

    p10 = cellTypeContr(trainedGP, "CD8", "Treg")
    p11 = cellTypeContr(trainedGP, "CD8", "Thelper")
    p12 = cellTypeContr(trainedGP, "CD8", "NK")
    p10rec = cellTypeContr(trainedGP, "CD8", "Treg", true)
    p11rec = cellTypeContr(trainedGP, "CD8", "Thelper", true)
    p12rec = cellTypeContr(trainedGP, "CD8", "NK", true)
    """p10CD4 = cellTypeContr(trainedGP4, "CD8", "Treg", false, true)
    p11CD4 = cellTypeContr(trainedGP4, "CD8", "Thelper", false, true)
    p12CD4 = cellTypeContr(trainedGP4, "CD8", "NK", false, true)
    p10CD4rec = cellTypeContr(trainedGP4, "CD8", "Treg", true, true)
    p11CD4rec = cellTypeContr(trainedGP4, "CD8", "Thelper", true, true)
    p12CD4rec = cellTypeContr(trainedGP4, "CD8", "NK", true, true)"""

    #draw(SVG("figureJ4.svg", 4000px, 1600px), p1)
    #draw(SVG("figureJ2.svg", 4000px, 1600px), gridstack([p13 p14 p15 p16; p17 p18 p19 p20; p21 p22 p23 p24]))
    #draw(SVG("figureJ4.svg", 4000px, 1600px), gridstack([p1 p1true]))

    #draw(SVG("figureJ4.svg", 2500px, 1600px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9; p10 p11 p12]))
    #draw(SVG("figureJ4.svg", 2500px, 6400px), gridstack([p1 p2 p3; p1CD4 p2CD4 p3CD4; p1rec p2rec p3rec; p1CD4rec p2CD4rec p3CD4rec; p4 p5 p6; p4CD4 p5CD4 p6CD4; p4rec p5rec p6rec; p4CD4rec p5CD4rec p6CD4rec; p7 p8 p9; p7CD4 p8CD4 p9CD4; p7rec p8rec p9rec; p7CD4rec p8CD4rec p9CD4rec; p10 p11 p12; p10CD4 p11CD4 p12CD4; p10rec p11rec p12rec; p10CD4rec p11CD4rec p12CD4rec]))

    draw(
        SVG("figureJ4.svg", 2500px, 3200px),
        gridstack([p1 p2 p3; p1rec p2rec p3rec; p4 p5 p6; p4rec p5rec p6rec; p7 p8 p9; p7rec p8rec p9rec; p10 p11 p12; p10rec p11rec p12rec]),
    )
end
