import LineSearches

dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")

""" Creates full vector of unknown values to be fit """
function getUnkVec()
    #kfwd, k4, k5, k16, k17, k22, k23, k27, endo, aendo, sort, krec, kdeg, k34, k35, k36, k37, k38, k39
    unkVecF = zeros(25)

    unkVecF[1] = 0.00125 # means of prior distributions from gc-cytokines paper
    unkVecF[2:7] .= 0.1
    unkVecF[8] = 1.0
    unkVecF[9] = 0.1
    unkVecF[10] = 0.678
    unkVecF[11] = 0.2
    unkVecF[12] = 0.01
    unkVecF[13] = 0.13
    unkVecF[14] = 2.0 # initial Treg stat
    unkVecF[15] = 0.5 # initial Thelp stat
    unkVecF[16] = 0.2 # initial NK stat
    unkVecF[17] = 0.2 # initial CD8 stat
    unkVecF[18:23] .= 0.001 # pSTAT Rates
    unkVecF[24] = 0.2
    unkVecF[25] = 0.2

    return unkVecF
end


""" Takes in full unkvec and constructs it into full fit parameters vector - TODO move this. """
function fitParams(ILs, unkVec::Vector{T}, recAbundances, CellType::String) where {T}
    kfbnd = 0.60
    paramvec = zeros(T, Nparams)
    paramvec[1:3] = ILs
    paramvec[4] = unkVec[1] #kfwd
    paramvec[5] = kfbnd * 10.0 #k1rev
    paramvec[6] = kfbnd * 144.0 #k2rev
    paramvec[7] = unkVec[2] #k4rev
    paramvec[8] = unkVec[3] #k5rev
    paramvec[9] = 12.0 * paramvec[8] / 1.5 #k10
    paramvec[10] = 63.0 * paramvec[8] / 1.5 #k11
    paramvec[11] = kfbnd * 0.065 #k13
    paramvec[12] = kfbnd * 438.0 #k14
    paramvec[13:16] = unkVec[4:7] #k15, k16, k17, k22
    paramvec[17] = kfbnd * 59.0 #k23
    paramvec[18] = unkVec[8] #k25
    paramvec[19] = 5.0 #endoadjust
    paramvec[20:24] = unkVec[9:13]
    paramvec[25:29] = receptor_expression(recAbundances, unkVec[9], unkVec[11], unkVec[12], unkVec[13])
    if CellType == "Treg"
        paramvec[30] = unkVec[14] #initial stat
    elseif CellType == "Thelper"
        paramvec[30] = unkVec[15] #initial stat
    elseif CellType == "NK"
        paramvec[30] = unkVec[16] #initial stat
    elseif CellType == "CD8"
        paramvec[30] = unkVec[17] #initial stat
    end
    paramvec[31:36] = unkVec[18:23]

    return paramvec
end


""" Adjusts the binding activity of ligand to match that of mutein """
function mutAffAdjust(paramVec::Vector{T}, dfRow) where {T}
    paramVec[5] = dfRow.IL2RaKD[1] * 0.6

    bgAdjust = (dfRow.IL2RBGKD[1] * 0.6) / paramVec[8]
    for ii in [6, 7, 8, 9, 10] # Adjust k2, k4, k5, k10, k11
        paramVec[ii] *= bgAdjust
    end

    return paramVec
end


""" Uses receptor abundance (from flow) and trafficking rates to calculate receptor expression rate at steady state. """
function receptor_expression(abundance, ke, kᵣ, sortF, kD)
    return abundance * ke / (1.0 + (kᵣ * (1.0 - sortF) / kD / sortF))
end


""" Calculates squared error for a given unkVec. """
function resids(x::Vector{T})::T where {T}
    @assert all(x .>= 0.0)
    df = importData()

    #get rid of IL15 and missing mutein
    df = df[df.Ligand .!= "R38Q/H16N", :]
    df = df[df.Ligand .!= "IL15", :]

    df.Time *= 60.0
    tps = unique(df.Time)
    sort!(tps)

    df.MeanPredict = similar(df.Mean, T)

    Threads.@threads for ligand in unique(df.Ligand)
        # Put the highest dose first so we catch a solving error early
        for dose in reverse(sort(unique(df.Dose)))
            ligVec = [dose, 0.0, 0.0]
            for cell in unique(df.Cell)
                idxx = findfirst(df.Cell .== cell)
                vector = vec(fitParams(ligVec, x, 10.0 .^ Vector{Float64}(df[idxx, [:IL15Ra, :IL2Ra, :IL2Rb, :IL7Ra, :gc]]), cell))
                if ligand != "IL2"
                    vector = mutAffAdjust(vector, df[findfirst(df.Ligand .== ligand), [:IL2RaKD, :IL2RBGKD]])
                end
                local yhat
                try
                    yhat = runCkine(tps, vector, pSTAT5 = true)
                catch e
                    if typeof(e) <: AssertionError
                        return Inf
                    else
                        rethrow(e)
                    end
                end

                for (ii, tt) in Iterators.enumerate(tps)
                    idxs = (df.Dose .== dose) .& (df.Time .== tt) .& (df.Ligand .== ligand) .& (df.Cell .== cell)
                    df[idxs, :MeanPredict] .= yhat[ii]
                end
            end
        end
    end

    #@assert all(df.MeanPredict .>= 0.0)
    dateFilt1 = filter(row -> string(row["Date"]) .== "2019-04-19", df)
    dateFilt2 = filter(row -> string(row["Date"]) .== "2019-05-02", df)
    dateFilt1.MeanPredict = dateFilt1.MeanPredict * x[24] * 1e6
    dateFilt2.MeanPredict = dateFilt2.MeanPredict * x[25] * 1e6
    df = vcat(dateFilt1, dateFilt2)
    #CSV.write("/home/brianoj/gcSolver.jl/data/fitTry.csv", x)
    # Convert relative scale.
    return norm((df.MeanPredict) - df.Mean)
end


""" Gets inital unkowns, optimizes them, and returns parameters of best fit"""
function runFit(; itern = 1000000)
    unk0 = invsoftplus.(getUnkVec())

    opts = Optim.Options(iterations = itern, show_trace = true)
    fit = optimize((x) -> resids(softplus.(x)), unk0, LBFGS(; linesearch = LineSearches.BackTracking()), opts, autodiff = :forward)

    @show fit

    return fit.minimizer
end

export getExpression, getUnkVec, fitParams
