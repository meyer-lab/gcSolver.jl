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
    #test = df.Ligand .== ligand
    #println("test = ", test)

    #println(first(df,5))
    
    #show(stdout, x)

    #first(realX, 5)
    realXtest = realX[1,:]
    println("realX = ", realXtest)
    
    #show(stdout, realXtest)

    realPreds = predict_f(gp, realX')

    biv = realX[1, 10]
    println("biv = ", biv)
    if biv == 1
        realX[:, 10] .= 0
    else
        realX[:, 10] .= 1
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

    if biv == 1
        predComp = gdf.plot(
            layer(predDF, x = :realPred, y = :fakePred, color = :cell, Geom.point),
            Guide.title(string("Ligand valency influence on (bivalent) ", ligand, " prediction")),
            Guide.xlabel("Correct Pred (log pSTAT, bivalent)"),
            Guide.ylabel("Incorrect Pred (log pSTAT, monovalent)"),
            Scale.color_discrete(),
            Guide.colorkey(title = "Cells"),
        )
    else
        predComp = gdf.plot(
            layer(predDF, x = :realPred, y = :fakePred, color = :cell, Geom.point),
            Guide.title(string("Ligand valency influence on (monovalent) ", ligand, " prediction")),
            Guide.xlabel("Correct Pred (log pSTAT, monovalent)"),
            Guide.ylabel("Incorrect Pred (log pSTAT, bivalent)"),
            Scale.color_discrete(),
            Guide.colorkey(title = "Cells"),
        )
    end
    return predComp
end

"""Generate plots"""
function figureJ6()
    x, y, df = gcSolver.getGPdata()
    trainedGP = gcSolver.gaussianProcess(x', y)

    p1 = BivContr(trainedGP, "IL2")
    p2 = BivContr(trainedGP, "IL15")
    p3 = BivContr(trainedGP, "WT N-term")
    p4 = BivContr(trainedGP, "H16N N-term")
    p5 = BivContr(trainedGP, "R38Q N-term")
    p6 = BivContr(trainedGP, "R38Q/H16N")
    p7 = BivContr(trainedGP, "V91K C-term")
    p8 = BivContr(trainedGP, "WT C-term")
    p9 = BivContr(trainedGP, "F42Q N-Term")

    #draw(SVG("figureJ6.svg", 3000px, 2000px), p1)
    draw(SVG("figureJ6.svg", 3000px, 2000px), gridstack([p1 p2 p3; p4 p5 p6; p7 p8 p9]))
end
