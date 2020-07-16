using GaussianProcesses
using Gadfly;
gdf = Gadfly;
import StatsBase: indicatormat


function getGPdata()
    fullData = importData()

    hotEnc = indicatormat(fullData.Cell)
    hotEncName = sort(unique(fullData.Cell))

    fullDataX = fullData[!, [:Dose, :Time, :IL2RaKD, :IL2RBGKD, :IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc]]
    fullDataY = log10.(fullData.Mean .+ 1.0)

    fullDataX[!, [:Dose, :IL2RaKD, :IL2RBGKD]] = log10.(fullDataX[!, [:Dose, :IL2RaKD, :IL2RBGKD]])

    for (ii, iiName) in enumerate(hotEncName)
        fullDataX[!, iiName] = vec(hotEnc[ii, :])
    end

    return Matrix{Float64}(fullDataX), vec(fullDataY), fullData
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
    CVplt = gdf.plot(
        layer(CVDF, x = :Yreal, y = :Y_pred, color = :Ligand, Geom.point),
        Guide.title(string("Leave-One-Mutein-Out CV")),
        Guide.xlabel("Actual pSTAT"),
        Guide.ylabel("Predicted pSTAT"),
        Scale.color_discrete(),
        Guide.colorkey(title = "Ligand"),
    )
    draw(SVG("LOMOcv.svg", 600px, 600px), CVplt)
    println(cor(y, y_pred))
end


function LOOcell()
    X, y, df = getGPdata()

    cells = unique(df.Cell)
    y_pred = zeros(length(y))
    cellList = Array{String}(undef, length(y))

    for cell in cells
        X_train = X[df.Cell .!== cell, :]
        y_train = y[df.Cell .!== cell]

        gp = gaussianProcess(X_train', y_train)

        yp, _ = predict_f(gp, X[df.Cell .== cell, :]')
        y_pred[df.Cell .== cell] .= yp
        cellList[df.Cell .== cell] .= cell
    end
    CVDF = DataFrame(Y_pred = y_pred, Yreal = y, Cell = cellList)

    CVplt = gdf.plot(
        layer(CVDF, x = :Yreal, y = :Y_pred, color = :Cell, Geom.point),
        Guide.title(string("Leave-One-Cell-Out CV")),
        Guide.xlabel("Actual pSTAT"),
        Guide.ylabel("Predicted pSTAT"),
        Scale.color_discrete(),
        Guide.colorkey(title = "Cell Type"),
    )
    draw(SVG("LOCOcv.svg", 600px, 600px), CVplt)
    println(cor(y, y_pred))
end

function cellHotEnc(cellType)
    if cellType == "CD8"
        return [1, 0, 0, 0]
    elseif cellType == "NK"
        return [0, 1, 0, 0]
    elseif cellType == "Thelper"
        return [0, 0, 1, 0]
    elseif cellType == "Treg"
        return [0, 0, 0, 1]
    else
        error("Unrecognized Cell Type")
    end
end

function runCkineVarPropGP(gp, xRow, sigma)::Vector

    # Sigma is the covariance matrix of the input parameters
    function jacF(x)

        pp = vcat(xRow[1:5], x[1:2], xRow[8], x[3], xRow[10:end])
        pp = reshape(pp, (13, 1))
        μ, σ² = predict_f(gp, pp)
        return μ
    end

    jac = zeros(3, 1)
    ForwardDiff.jacobian!(jac, jacF, append!(xRow[6:7], xRow[9]))

    # Just return the diagonal for the marginal variance
    return diag(transpose(jac) * sigma * jac)
end
