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


function figureJ1()
    print(responseDF[1:10])
    p1 = trialplot()
    p2 = trialplot()
    p3 = trialplot()
    p4 = trialplot()

    draw(SVG("figureJ1.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end
