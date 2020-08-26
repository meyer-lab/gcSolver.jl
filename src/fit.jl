import LineSearches: BackTracking

dataDir = joinpath(dirname(pathof(gcSolver)), "..", "data")

""" Creates full vector of unknown values to be fit """
function getUnkVec()
    #kfwd, k4, k5, k16, k17, k22, k23, k27, endo, aendo, sort, krec, kdeg, k34, k35, k36, k37, k38, k39
    p = 0.1ones(23)

    p[1] = 0.001 # means of prior distributions from gc-cytokines paper
    p[8] = 1.0
    p[10] = 0.678
    p[11] = 0.2
    p[12] = 0.01
    p[13] = 0.13
    p[14] = 2.0 # initial Treg stat
    p[15] = 0.5 # initial Thelp stat
    p[16] = 0.2 # initial NK stat
    p[17] = 0.2 # initial CD8 stat
    p[18:23] .= 0.001 # pSTAT Rates

    return p
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
function mutAffAdjust(p::Vector{T}, dfRow) where {T}
    p[5] = dfRow.IL2RaKD[1] * 0.6
    p[6:10] .*= p[5] / p[8]
    return p
end


""" Uses receptor abundance (from flow) and trafficking rates to calculate receptor expression rate at steady state. """
function receptor_expression(abundance, ke, kᵣ, sortF, kD)
    return abundance * ke / (1.0 + (kᵣ * (1.0 - sortF) / kD / sortF))
end


""" Calculates squared error for a given unkVec. """
function resids(x::Vector{T})::T where {T}
    @assert all(x .>= 0.0)
    df = importData()

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

    # @assert all(df.MeanPredict .>= 0.0)
    dateFilt1 = filter(row -> string(row["Date"]) .== "4/19/2019", df)
    dateFilt2 = filter(row -> string(row["Date"]) .== "5/2/2019", df)
    # Convert relative scale.
    dateFilt1.MeanPredict .*= dateFilt1.MeanPredict \ dateFilt1.Mean
    dateFilt2.MeanPredict .*= dateFilt2.MeanPredict \ dateFilt2.Mean
    return norm(dateFilt1.MeanPredict - dateFilt1.Mean) + norm(dateFilt2.MeanPredict - dateFilt2.Mean)
end


""" Gets inital unkowns, optimizes them, and returns parameters of best fit"""
function runFit(; itern = 1000000)
    x₀ = invsoftplus.(getUnkVec())

    opts = Optim.Options(iterations = itern, show_trace = true)
    fit = optimize((x) -> resids(softplus.(x)), x₀, LBFGS(; linesearch = BackTracking()), opts, autodiff = :forward)

    @show fit

    return fit.minimizer
end

export getExpression, getUnkVec, fitParams
