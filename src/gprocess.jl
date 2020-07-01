using GaussianProcesses
import StatsBase: indicatormat


function getGPdata()
    yData = getyVec()
    affDF = DataFrame!(CSV.File(joinpath(dataDir, "mutAffData.csv")))

    fullData = innerjoin(yData, affDF, on = :Ligand => :Mutein)

    # Hot encode cell types
    hotEnc = indicatormat(fullData.Cell)
    hotEncName = sort(unique(fullData.Cell))

    fullDataX = fullData[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD]]
    fullDataY = log10.(fullData.Mean .+ 1.0)

    fullDataX[!, :Dose] = log10.(fullDataX[!, :Dose])
    fullDataX[!, :IL2RaKD] = log10.(fullDataX[!, :IL2RaKD])
    fullDataX[!, :IL2RBGKD] = log10.(fullDataX[!, :IL2RBGKD])

    # Encode the cell types with one-hot
    for (ii, iiName) in enumerate(hotEncName)
        fullDataX[!, iiName] = vec(hotEnc[ii, :])
    end

    return Matrix(fullDataX), vec(fullDataY), fullData
end


" Assemble Gaussian process model. "
function gaussianProcess(X, y::Vector)
    mZero = MeanZero()
    kern = SE(0.0, 0.0)

    gp = GP(X, y, mZero, kern)

    optimize!(gp)

    return gp
end


function gaussianTest()
    X, y, _ = getGPdata()

    gp = gaussianProcess(X', y)

    μ, σ2 = predict_f(gp, X[1:2, :]')
end


function LOOmutein()
    X, y, df = getGPdata()

    muteins = unique(df.Ligand)
    y_pred = zeros(length(y))

    for mutein in muteins
        X_train = X[df.Ligand .!== mutein, :]
        y_train = y[df.Ligand .!== mutein]

        gp = gaussianProcess(X_train', y_train)

        yp, _ = predict_f(gp, X[df.Ligand .== mutein, :]')
        y_pred[df.Ligand .== mutein] .= yp
    end

    println(cor(y, y_pred))
end
