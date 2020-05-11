using Flux
using Flux.Data: DataLoader
using Flux: @epochs
using DataFrames
import StatsBase: indicatormat
import NNlib


function build_model()
    return Chain(Dense(8, 256, NNlib.sigmoid), Dropout(0.05), Dense(256, 1))
end


function getMLPdata()
    yData = getyVec()
    affDF = CSV.read(joinpath(dataDir, "mutAffData.csv"), copycols = true)

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

    return Float32.(transpose(Matrix(fullDataX))), Float32.(fullDataY)
end


function loss_all(XX, yy, model)
    println(cor(vec(model(XX)), yy))

    return norm(model(XX) .- yy)
end


function train()
    # Load Data
    XX, y = getMLPdata()

    train_loader = DataLoader(XX, y, batchsize = length(y) / 2, shuffle = true)

    # Construct model
    m = build_model()
    loss(xIn, yyIn) = norm(m(xIn) .- yyIn)
    evalcb = () -> @show(loss_all(XX, y, m))

    ## Training
    # Keep getting ~1264
    @epochs 40000 Flux.train!(loss, params(m), train_loader, ADAM(), cb = evalcb)

    return y, m(XX)
end
