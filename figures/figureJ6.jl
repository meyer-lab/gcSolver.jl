""" This file builds the depletion manuscript, Figure 4. """

using Gadfly;
using gcSolver;
using DataFrames;
using GaussianProcesses;
gdf = Gadfly;

function BivContr(gp, ligand)
    x, y, df = gcSolver.getGPdata()

    predDF = DataFrame(realPred = Float64[], fakePred = Float64[], cell = String[])

    realX = x[df.Ligand .== ligand, :]

    realPreds = predict_f(gp, realX')
    biv = realX[1, size(realX, 2)]
    if biv == 1
        realX[:, size(realX, 2)] .= 0
    else
        realX[:, size(realX, 2)] .= 1
    end

    compPreds = predict_f(gp, realX')

    cells = df[df.Ligand .== ligand, :].Cell

    # first tuple returned by predict_f is the predictions and the second tuple returned is the standard deviations
    realVals = realPreds[1]
    fakeVals = compPreds[1]
    #for i in enumerate(realVals)
    for i = 1:length(cells)
        real = realVals[i]
        fake = fakeVals[i]
        cel = cells[i]
        push!(predDF, (real, fake, cel))
    end


    predComp = gdf.plot(
        layer(predDF, x = :realPred, y = :fakePred, color = :cell, Geom.point),
        Guide.title(string("Ligand valency influence on ", ligand, " prediction")),
        Guide.xlabel("Correct Pred"),
        Guide.ylabel("Incorrect Pred"),
        Scale.color_discrete(),
        Guide.colorkey(title = "Cells"),
    )
    return predComp
end

"""Generate plots"""
function figureJ6()
    x, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(x', y)

    #x4, y4, df4 = getGPdata(true)
    #trainedGP4 = gaussianProcess(x4', y4)

    p1 = BivContr(trainedGP, "IL2")
    p2 = BivContr(trainedGP, "IL15")
    p3 = BivContr(trainedGP, "WT N-term")
    p4 = BivContr(trainedGP, "H16N N-term")
    p5 = BivContr(trainedGP, "R38Q N-term")
    p6 = BivContr(trainedGP, "R38Q/H16N")

    draw(SVG("figureJ6.svg", 3000px, 2000px), gridstack([p1 p2 p3; p4 p5 p6]))
end
