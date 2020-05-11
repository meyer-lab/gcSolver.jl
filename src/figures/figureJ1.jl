""" This file builds the depletion manuscript, Figure 1. """

""" Plot an example isobologram. """
function trialplot()
    X = [1, 2, 3]
    Y = [1, 2, 3]
    pl = plot(x=X, y=Y);
    return pl
end


function figureJ1()
    p1 = trialplot()
    p2 = trialplot()
    p3 = trialplot()
    p4 = trialplot()

    draw(SVG("figureJ1.svg", 1000px, 800px), gridstack([p1 p2; p3 p4]))
end
