using Bijectors
using Optim
using Statistics


function optimizefunc(flParams)
    #Fits your flow to a target distribution (baseDist). This gives the opimized flow to transform your base distribution to a multivariate normal distribution.
    baseDist = MvNormal([0, 0], ones(2))
    fitflow = PlanarLayer([flParams[1], flParams[2]], [flParams[3], flParams[4]], [flParams[5]]) #∘ RadialLayer([x[6]], [x[7]], [x[8], x[9]])
    flofitT = transformed(baseDist, fitflow)
    answer = -sum(logpdf(flofitT, targetsamps))
    return answer
end


function getpoints(fitMin, npoints)
    #Uses inverse of fit flow to transform a distribution to match your original points in space.
    baseDist = MvNormal(zeros(2), ones(2))
    fitflow = PlanarLayer([fitMin[1], fitMin[2]], [fitMin[3], fitMin[4]], [fitMin[5]]) #∘ RadialLayer([x[6]], [x[7]], [x[8], x[9]])
    fitflow = inv(fitflow)
    answer = rand(baseDist, npoints)
    answer = fitflow(answer)
    return answer
end


function fitFlow()
    #uses optimizefunc to return optimized flow parameters.
    result = optimize(optimizefunc, rand(5), LBFGS(), options)
    cost = result.minimum
    fitParams = result.minimizer
    return fitParams, cost
end
