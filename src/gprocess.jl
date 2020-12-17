using GaussianProcesses, StatsBase

function getGPdata(log=true)
    fullData = importData()

    hotEnc = indicatormat(fullData.Cell)
    hotEncName = sort(unique(fullData.Cell))

    fullDataX = fullData[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc, :Bivalent]]
    fullDataY = log10.(fullData.Mean .+ 1.0)
    if log
        fullDataY = log10.(fullData.Mean .+ 1.0)
    else
        fullDataY = fullData.Mean
    end

    fullDataX[!, [:Dose, :IL2RaKD, :IL2RBGKD]] = log10.(fullDataX[!, [:Dose, :IL2RaKD, :IL2RBGKD]])

    for (ii, iiName) in enumerate(hotEncName)
        fullDataX[!, iiName] = vec(hotEnc[ii, :])
    end

    combCD4 = fullDataX[:, size(fullDataX, 2)]
    for iii = 1:length(combCD4)
        if combCD4[iii] == 1
            fullDataX[iii, 13] = 1
        end
    end

    Xdat = DataFrame(Matrix{Float64}(fullDataX))
    Ydat = DataFrame(Y= fullDataY)
    fullDataX = fullDataX[:, 1:13]

    return Matrix{Float64}(fullDataX), vec(fullDataY), fullData
end


" Assemble Gaussian process model. "
function gaussianProcess(X, y::Vector)
    mZero = MeanZero()

    lscales = zeros(Float64, size(X)[1])
    kern = RQ(lscales, 0.0, 0.0)
    indices = collect(1:size(X)[2])
    idxs = sample(indices, convert(Int64, round(size(X)[2]/5, digits=0)), replace=false)
    Xtrain = X[:, idxs]
    Ytrain = y[idxs]
    gp = GPE(Xtrain, Ytrain, mZero, kern)

    optimize!(gp)

    return gp
end


function gaussianTest()
    X, y, _ = getGPdata()

    gp = gaussianProcess(X', y)

    μ, σ2 = predict_f(gp, X[1:2, :]')
end



function LOOCV()
    X, y, df = getGPdata()

    muteins = unique(df.Ligand)
    y_pred = zeros(length(y))
    muteinList = Array{String}(undef, length(y))
    gp = gaussianProcess(X', y)

    y_pred, _ = GaussianProcesses.predict_LOO(gp)
    muteinList = df.Ligand
    CVDF = DataFrame(Y_pred = y_pred, Yreal = y, Ligand = muteinList)

    println(cor(y, y_pred))
end


function LOOmutein()
    X, y, df = getGPdata()

    muteins = unique(df.Ligand)
    y_pred = zeros(length(y))
    muteinList = Array{String}(undef, length(y))

    for mutein in muteins
        X_train = X[df.Ligand .!== mutein, :]
        y_train = y[df.Ligand .!== mutein]

        gp = gaussianProcess(X_train', y_train)

        yp, _ = predict_f(gp, X[df.Ligand .== mutein, :]')
        y_pred[df.Ligand .== mutein] .= yp
        muteinList[df.Ligand .== mutein] .= mutein
    end
    CVDF = DataFrame(Y_pred = y_pred, Yreal = y, Ligand = muteinList)

    println(cor(y, y_pred))
end


function LOOcell()
    X, y, df = getGPdata()

    y_pred = zeros(length(y))
    cellList = Array{String}(undef, length(y))

    for cell in unique(df.Cell)
        X_train = X[df.Cell .!== cell, :]
        y_train = y[df.Cell .!== cell]

        gp = gaussianProcess(X_train', y_train)

        yp, _ = predict_f(gp, X[df.Cell .== cell, :]')
        y_pred[df.Cell .== cell] .= yp
        cellList[df.Cell .== cell] .= cell
    end
    CVDF = DataFrame(Y_pred = y_pred, Yreal = y, Cell = cellList)

    println(cor(y, y_pred))
end

function cellHotEnc(cellType)
    if cellType == "CD8"
        return [1, 0, 0]
    elseif cellType == "NK"
        return [0, 1, 0]
    elseif cellType == "Thelper"
        return [0, 0, 1]
    elseif cellType == "Treg"
        return [0, 0, 1]
    else
        error("Unrecognized Cell Type")
    end
end

function runCkineVarPropGP(gp, xRow, sigma, cov = false)::Vector

    # Sigma is the covariance matrix of the input parameters
    if cov
        #take only variance in pstat explained by IL2Ra variance
        function jacFCov(x)
            pp = vcat(xRow[1:5], x[1], xRow[7:end])
            pp = reshape(pp, size(xRow, 1), 1)
            μ, σ² = predict_f(gp, pp)
            return μ
        end

        jac = zeros(1, 1)
        ForwardDiff.jacobian!(jac, jacFCov, [xRow[6]])

        # Just return the diagonal for the marginal variance
        return diag(transpose(jac) * sigma[1, 1] * jac)

    else
        #variance explained by all receptor expression discrepencies
        function jacF(x)

            pp = vcat(xRow[1:5], x[1:2], xRow[8], x[3], xRow[10:end])
            pp = reshape(pp, size(xRow, 1), 1)
            μ, σ² = predict_f(gp, pp)
            return μ
        end

        jac = zeros(3, 1)
        ForwardDiff.jacobian!(jac, jacF, append!(xRow[6:7], xRow[9]))

        # Just return the diagonal for the marginal variance
        return diag(transpose(jac) * sigma * jac)
    end
end
