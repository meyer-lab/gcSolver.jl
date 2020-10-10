using Distributed
import LineSearches: InitialStatic
using CSV

dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")

""" Creates full vector of unknown values to be fit """
function getUnkVec()
    #kfwd, k4, k5, k16, k17, k22, k23, k27, endo, aendo, sort, krec, kdeg, k34, k35, k36, k37, k38, k39
    p = 0.1ones(23)

    p[1] = 0.001 # means of prior distributions from gc-cytokines paper
    p[8] = 1.0 #27
    p[10] = 0.678
    p[11] = 0.2
    p[12] = 0.01
    p[13] = 0.13
    p[14] = 5.0 # initial Treg stat
    p[15] = 0.2 # initial Thelp stat
    p[16] = 0.3 # initial NK stat
    p[17] = 0.1 # initial CD8 stat
    p[18:23] .= 0.001

    return p
end


""" Takes in full unkvec and constructs it into full fit parameters vector - TODO move this. """
function fitParams(ILs, unkVec::Vector{T}, recAbundances, CellType::String) where {T}
    kfbnd = 0.60
    farhatVec = getFarhatVec()
    paramvec = zeros(T, Nparams)
    paramvec[1:3] = ILs
    paramvec[4] = farhatVec[1] #unkVec[1] #kfwd
    paramvec[5] = kfbnd * 10.0 #k1rev
    paramvec[6] = kfbnd * 144.0 #k2rev
    paramvec[7] = farhatVec[2] # unkVec[2] #k4rev
    paramvec[8] = farhatVec[3]#unkVec[3] #k5rev
    paramvec[9] = 12.0 * paramvec[8] / 1.5 #k10
    paramvec[10] = 63.0 * paramvec[8] / 1.5 #k11
    paramvec[11] = kfbnd * 0.065 #k13
    paramvec[12] = kfbnd * 438.0 #k14
    paramvec[13:16] = farhatVec[4:7]# unkVec[4:7] #k16, k17, k22, k23
    paramvec[17] = kfbnd * 59.0 #k25
    paramvec[18] = farhatVec[8]#unkVec[8]  #k27
    paramvec[19] = 5.0 #endoadjust #make sure this is right
    ke = farhatVec[9]#unkVec[9] 
    aendo = farhatVec[10]#unkVec[10]
    sortF = farhatVec[11]# tanh(unkVec[11]) * 0.95 + 0.01
    krec = farhatVec[12]#unkVec[12]
    kdeg = farhatVec[13] #unkVec[13]
    paramvec[20] = ke
    paramvec[21] = aendo
    paramvec[22] = sortF
    paramvec[23] = krec
    paramvec[24] = kdeg
    paramvec[25:29] = recAbundances * ke / (1.0 + (krec * (1.0 - sortF) / kdeg / sortF))
    if CellType == "Treg"
        paramvec[30] = unkVec[1] #initial stat
    elseif CellType == "Thelper"
        paramvec[30] = unkVec[2] #initial stat
    elseif CellType == "NK"
        paramvec[30] = unkVec[3] #initial stat
    elseif CellType == "CD8"
        paramvec[30] = unkVec[4] #initial stat
    else
        @assert false
    end
    paramvec[31:36] = unkVec[5:10]

    return paramvec
end


""" Adjusts the binding activity of ligand to match that of mutein """
function mutAffAdjust(p::Vector{T}, dfRow) where {T}
    p[5] = dfRow.IL2RaKD[1] * 0.6
    bgAdjust = (dfRow.IL2RBGKD[1] * 0.6) / p[8]
    p[6:10] .*= bgAdjust
    return p
end


""" Calculates squared error for a given unkVec. """
function resids(x::Vector{T})::T where {T}
    @assert all(x .>= 0.0)
    df = importData(true)
    df = df[df.Ligand .!= "IL15", :]

    sort!(df, :Time)
    df.Time *= 60.0

    df.MeanPredict = similar(df.Mean, T)
    FutureDict = Dict()
    cost = 0.0

    for ligand in unique(df.Ligand)
        # Put the highest dose first so we catch a solving error early
        for dose in reverse(sort(unique(df.Dose)))
            if ligand == "IL15"
                ligVec = [0.0, dose, 0.0]
            else
                ligVec = [dose, 0.0, 0.0]
            end

            for cell in unique(df.Cell)
                idxx = findfirst(df.Cell .== cell)
                recpE = 10.0 .^ Vector{Float64}(df[idxx, [:IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc]])
                vector = vec(fitParams(ligVec, x, recpE, cell))
                if ligand !== "IL15" && ligand !== "IL2"
                    vector = mutAffAdjust(vector, df[findfirst(df.Ligand .== ligand), [:IL2RaKD, :IL2RBGKD]])
                end
                idxs = (df.Dose .== dose) .& (df.Ligand .== ligand) .& (df.Cell .== cell)

                tpss = df[idxs, :Time]
                # Make sure duplicate times are not considered duplicates
                tpss += range(0.0, 0.01; length = length(tpss))

                # Regularize for exploding values
                cost += sum(softplus.(vector .- 1.0e6))

                FutureDict[(dose, ligand, cell)] = @spawnat :any runCkine(tpss, vector; pSTAT5 = true)
            end
        end
    end

    for (key, val) in FutureDict
        idxs = (df.Dose .== key[1]) .& (df.Ligand .== key[2]) .& (df.Cell .== key[3])

        df[idxs, :MeanPredict] = fetch(val)
    end

    @assert all(df.MeanPredict .>= -0.01)

    for date in unique(df.Date)
        dateFilt = filter(row -> string(row["Date"]) .== date, df)
        dateFilt.MeanPredict .*= dateFilt.MeanPredict \ dateFilt.Mean
        cost += norm(dateFilt.MeanPredict - dateFilt.Mean)
    end

    return cost
end


""" Gets inital unkowns, optimizes them, and returns parameters of best fit"""
function runFit(; itern = 1000000)
    x₀ = invsoftplus.(getUnkVec()[14:23])

    opts = Optim.Options(iterations = itern, show_trace = true, extended_trace = true, allow_f_increases = true)
    lsi = InitialStatic(; alpha = 0.0001)
    fit = optimize((x) -> resids(softplus.(x)), x₀, LBFGS(; m = 100, alphaguess = lsi), opts, autodiff = :forward)

    @show fit

    return fit.minimizer
end

export getExpression, getUnkVec, fitParams, runCkine


""" Writes a CSV for conversion factors between predictions and pSTAT values """
function getDateConvDict()
    fitVec = importFit()
    fitVec = convert(Vector{Float64}, fitVec[!, :Fit])
    x = softplus.(fitVec)

    df = importData(true)

    sort!(df, :Time)
    df.Time *= 60.0

    df.MeanPredict = similar(df.Mean, T)
    FutureDict = Dict()

    for ligand in unique(df.Ligand)
        # Put the highest dose first so we catch a solving error early
        for dose in reverse(sort(unique(df.Dose)))
            if ligand == "IL15"
                ligVec = [0.0, dose, 0.0]
            else
                ligVec = [dose, 0.0, 0.0]
            end

            for cell in unique(df.Cell)
                idxx = findfirst(df.Cell .== cell)
                recpE = 10.0 .^ Vector{Float64}(df[idxx, [:IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc]])
                vector = vec(fitParams(ligVec, x, recpE, cell))

                if ligand != "IL15"
                    vector = mutAffAdjust(vector, df[findfirst(df.Ligand .== ligand), [:IL2RaKD, :IL2RBGKD]])
                end
                idxs = (df.Dose .== dose) .& (df.Ligand .== ligand) .& (df.Cell .== cell)

                tpss = df[idxs, :Time]
                # Make sure duplicate times are not considered duplicates
                tpss += range(0.0, 0.01; length = length(tpss))

                # Regularize for exploding values

                FutureDict[(dose, ligand, cell)] = @spawnat :any runCkine(tpss, vector; pSTAT5 = true)
            end
        end
    end

    for (key, val) in FutureDict
        idxs = (df.Dose .== key[1]) .& (df.Ligand .== key[2]) .& (df.Cell .== key[3])

        df[idxs, :MeanPredict] = fetch(val)
    end

    @assert all(df.MeanPredict .>= -0.01)

    dateConvDF = DataFrame(Date = [], Conv = [])
    for date in unique(df.Date)
        dateFilt = filter(row -> string(row["Date"]) .== date, df)
        convFact = dateFilt.MeanPredict \ dateFilt.Mean
        miniDF = DataFrame(Date = date, Conv = convFact)
        dateFilt.MeanPredict .*= dateFilt.MeanPredict \ dateFilt.Mean
        append!(dateConvDF, miniDF)
    end


    CSV.write(joinpath(dataDir, "DateConvFrame.csv"), dateConvDF)
end


""" Writes results of Farhat et al. model fitting to unknown vector """
function getFarhatVec()
    #kfwd, k4, k5, k16, k17, k22, k23, k27, endo, aendo, sort, krec, kdeg, k34, k35, k36, k37, k38, k39
    p = 0.1ones(23)
    p[1] = 0.003 #0.001 # means of prior distributions from gc-cytokines paper
    p[2] = 12.25
    p[3] = 0.08
    p[4] = 28
    p[5] = 0.04
    p[6] = 4.57
    p[7] = 25.12
    p[8] = 1.0
    p[9] = 0.09
    p[10] = 2.5
    p[11] = 0.125
    p[12] = 0.009
    p[13] = 0.0068

    p[14] = 2.0 # initial Treg stat
    p[15] = 0.5 # initial Thelp stat
    p[16] = 0.2 # initial NK stat
    p[17] = 0.2 # initial CD8 stat
    p[18:23] .= 0.001 # pSTAT Rates
    return p
end
