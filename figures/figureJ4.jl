""" This file builds the depletion manuscript, Figure 4. """

using Gadfly;
using gcSolver;
using DataFrames;
using GaussianProcesses;
gdf = Gadfly;

function cellTypeContr(gp, realType, compType, recExp = false, biv = true)

    x, y, df = gcSolver.getGPdata()

    predDF = DataFrame(realPred = Float64[], fakePred = Float64[], ligand = String[])

    #RealX = only rows with x axis cell-type, want it to make predictions
    realX = x[df.Cell .== realType, :]

    realPreds = predict_f(gp, realX')
    # returns tuple of arrays

    if recExp == true
        compExpression = x[df.Cell .== compType, :][1, 5:9]
        realX[:, 5:9] = repeat(compExpression , outer = [1, size(realX)[1]])'
    else
        hotEnc = gcSolver.cellHotEnc(compType)
        realX[:, 11:size(realX, 2)] = repeat(hotEnc, outer = [1, size(realX)[1]])'
    end

    compPreds = predict_f(gp, realX')

    ligands = df[df.Cell .== realType, :]
    if biv == true
        ligands = filter(row -> row["Bivalent"] .== biv, ligands)
    end
    ligands = ligands.Ligand
    # first tuple returned by predict_f is the predictions and the second tuple returned is the standard deviations
    realVals = realPreds[1]
    fakeVals = compPreds[1]

    for i = 1:length(ligands)
        real = realVals[i]
        fake = fakeVals[i]
        lig = ligands[i]
        push!(predDF, (real, fake, lig))
    end

    if biv == true
        if recExp == true
            predComp = gdf.plot(
                layer(predDF, x = :realPred, y = :fakePred, color = :ligand, Geom.point),
                Guide.title(string(compType, " influence on ", realType, " prediction (Receptor profile, CD4+)")),
                Guide.xlabel("Correct Pred"),
                Guide.ylabel("Incorrect Pred"),
                Scale.color_discrete(),
                Guide.colorkey(title = "All Muteins"),
            )
        else
            predComp = gdf.plot(
                layer(predDF, x = :realPred, y = :fakePred, color = :ligand, Geom.point),
                Guide.title(string(compType, " influence on ", realType, " prediction (CD4+)")),
                Guide.xlabel("Correct Pred"),
                Guide.ylabel("Incorrect Pred"),
                Scale.color_discrete(),
                Guide.colorkey(title = "All Muteins"),
            )
        end
    else
        if recExp == true
            predComp = gdf.plot(
                layer(predDF, x = :realPred, y = :fakePred, color = :ligand, Geom.point),
                Guide.title(string(compType, " influence on ", realType, " prediction (Receptor profile, CD4+)")),
                Guide.xlabel("Correct Pred"),
                Guide.ylabel("Incorrect Pred"),
                Scale.color_discrete(),
                Guide.colorkey(title = "Monomeric Muteins"),
            )
        else
            predComp = gdf.plot(
                layer(predDF, x = :realPred, y = :fakePred, color = :ligand, Geom.point),
                Guide.title(string(compType, " influence on ", realType, " prediction (CD4+)")),
                Guide.xlabel("Correct Pred"),
                Guide.ylabel("Incorrect Pred"),
                Scale.color_discrete(),
                Guide.colorkey(title = "Monomeric Muteins"),
            )
        end
    end


    return predComp
end

"""Generate plots"""
function figureJ4()
    x, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(x', y)

    """p1 = cellTypeContr(trainedGP, "Treg", "Thelper", false)
    p2 = cellTypeContr(trainedGP, "Treg", "NK", false)
    p3 = cellTypeContr(trainedGP, "Treg", "CD8", false)
    p1rec = cellTypeContr(trainedGP, "Treg", "Thelper", true)
    p2rec = cellTypeContr(trainedGP, "Treg", "NK", true)
    p3rec = cellTypeContr(trainedGP, "Treg", "CD8", true)

    p4 = cellTypeContr(trainedGP, "Thelper", "Treg")
    p5 = cellTypeContr(trainedGP, "Thelper", "NK")
    p6 = cellTypeContr(trainedGP, "Thelper", "CD8")
    p4rec = cellTypeContr(trainedGP, "Thelper", "Treg", true)
    p5rec = cellTypeContr(trainedGP, "Thelper", "NK", true)
    p6rec = cellTypeContr(trainedGP, "Thelper", "CD8", true)
    
    p7 = cellTypeContr(trainedGP, "NK", "Treg")
    p8 = cellTypeContr(trainedGP, "NK", "Thelper")
    p9 = cellTypeContr(trainedGP, "NK", "CD8")
    p7rec = cellTypeContr(trainedGP, "NK", "Treg", true)
    p8rec = cellTypeContr(trainedGP, "NK", "Thelper", true)
    p9rec = cellTypeContr(trainedGP, "NK", "CD8", true)

    p10 = cellTypeContr(trainedGP, "CD8", "Treg")
    p11 = cellTypeContr(trainedGP, "CD8", "Thelper")
    p12 = cellTypeContr(trainedGP, "CD8", "NK")
    p10rec = cellTypeContr(trainedGP, "CD8", "Treg", true)
    p11rec = cellTypeContr(trainedGP, "CD8", "Thelper", true)
    p12rec = cellTypeContr(trainedGP, "CD8", "NK", true)"""

    p1 = cellTypeContr(trainedGP, "Treg", "Thelper", false, 0)
    p2 = cellTypeContr(trainedGP, "Treg", "Thelper", false, 1)
    p3 = cellTypeContr(trainedGP, "Treg", "Thelper", true, 0)
    p4 = cellTypeContr(trainedGP, "Treg", "Thelper", true, 1)

    draw(SVG("figureJ4.svg", 4000px, 1600px), gridstack([p1 p2; p3 p4]))
    
    """draw(
        SVG("figureJ4.svg", 2500px, 3200px),
        gridstack([p1 p2 p3; p1rec p2rec p3rec; p4 p5 p6; p4rec p5rec p6rec; p7 p8 p9; p7rec p8rec p9rec; p10 p11 p12; p10rec p11rec p12rec]),
    )"""
end
