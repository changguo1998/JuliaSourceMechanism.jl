function pick_windowratio(x::AbstractVector{<:Real}, wl::Integer)
    L = length(x)
    r = zeros(L - 2 * wl)
    for i in eachindex(r)
        r[i] = std(x[i+wl:i+2*wl]) / std(x[i:i+wl])
    end
    (_, j) = findmax(r)
    return (j + wl, r)
end

_procfunc(x) = mean(abs, x)

function _stalta(x::AbstractVector{<:Real}, ws::Integer, wl::Integer)
    if ws > wl
        error("short window larger than long window")
    end
    L = length(x)
    r = zeros(L-2*wl)
    for i = eachindex(r)
        r[i] = _procfunc(x[i+wl:i+wl+ws])/_procfunc(x[i:i+wl])
    end
    return r
end

function pick_stalta(x::AbstractVector{<:Real}, ws::Integer, wl::Integer)
    r = _stalta(x, ws, wl)
    return (argmax(r)+wl, maximum(r))
end

function _freedom(x::AbstractVecOrMat)
    if size(x, 2) == 1
        return 1.0
    else
        F = svd(x)
        return sum(F.S) / maximum(F.S)
    end
end

function pick_freedom(x::AbstractVecOrMat{<:Real}, wl::Integer)
    L = size(x, 1)
    r = zeros(L - wl)
    for i in eachindex(r)
        r[i] = _freedom(x[i:i+wl, :])
    end
    (_, j) = findmin(r)
    return (j + wl, r)
end

function detrendandtaper!(x::AbstractVecOrMat)
    for v in eachcol(x)
        SeisTools.DataProcess.detrend!(v)
        SeisTools.DataProcess.taper!(v)
    end
    return nothing
end

"""
dc2ts(sdr::Vector) where T <: Real -> m::Vector{AbstractFloat}

convert double couple source to moment tensor in NED coordinary
"""
function dc2ts(sdr::Vector{T}) where {T<:Real}
    s = deg2rad(sdr[1])
    d = deg2rad(sdr[2])
    r = deg2rad(sdr[3])
    m = zeros(6)
    m[1] = -1 * (sin(2 * s) * sin(d) * cos(r) + (sin(s))^2 * sin(2 * d) * sin(r))
    m[2] = sin(2 * s) * sin(d) * cos(r) - (cos(s))^2 * sin(2 * d) * sin(r)
    m[3] = sin(2 * d) * sin(r)
    m[4] = cos(2 * s) * sin(d) * cos(r) + 0.5 * sin(2 * s) * sin(2 * d) * sin(r)
    m[5] = -1 * (cos(s) * cos(d) * cos(r) + sin(s) * cos(2 * d) * sin(r))
    m[6] = -1 * (sin(s) * cos(d) * cos(r) - cos(s) * cos(2 * d) * sin(r))
    return m
end

function calcgreen!(env::Setting; showinfo::Bool=false)
    taglist = String[]
    idxlist = Int[]
    for i in eachindex(env["stations"])
        s = env["stations"][i]
        tag = String(s["network"] * "." * s["station"])
        if !(tag in taglist)
            push!(taglist, tag)
            push!(idxlist, i)
        end
    end
    Threads.@threads for i in idxlist
        s = env["stations"][i]
        (dist, az, _) = SeisTools.Geodesy.distance(env["event"]["latitude"], env["event"]["longitude"], s["meta_lat"],
                                                   s["meta_lon"])
        s["base_distance"] = dist
        s["base_azimuth"] = az
        Green.calc(s, env; showinfo=showinfo)
    end
    return nothing
end

function loaddata!(env::Setting; showinfo::Bool=false)
    for s in env["stations"]
        (dist, az, _) = SeisTools.Geodesy.distance(env["event"]["latitude"], env["event"]["longitude"], s["meta_lat"],
                                                   s["meta_lon"])
        s["base_distance"] = dist
        s["base_azimuth"] = az
    end
    for s in env["stations"]
        t = SeisTools.SAC.read(normpath(env["dataroot"], "sac", s["meta_file"]))
        trim_bt = s["base_trim"][1]
        trim_et = s["base_trim"][2]
        (sbt, tw, _) = SeisTools.DataProcess.cut(t.data, SeisTools.SAC.DateTime(t.hdr), trim_bt, trim_et,
                                                 Millisecond(round(Int, t.hdr["delta"] * 1000)))
        s["base_begintime"] = sbt
        s["base_record"] = tw
        Green.load!(s, env; showinfo=showinfo)
        detrendandtaper!(s["base_record"])
        detrendandtaper!(s["green_fun"])
    end
    return nothing
end

"""
preprocess!(env::Setting, modules::Vector{Module}; showinfo::Bool=false)
"""
function preprocess!(env::Setting, modules::Vector{Module}; showinfo::Bool=false)
    for s in env["stations"]
        t = SeisTools.SAC.read(normpath(env["dataroot"], "sac", s["meta_file"]))
        trim_bt = s["base_trim"][1]
        trim_et = s["base_trim"][2]
        (sbt, tw, _) = SeisTools.DataProcess.cut(t.data, SeisTools.SAC.DateTime(t.hdr), trim_bt, trim_et,
                                                 Millisecond(round(Int, t.hdr["delta"] * 1000)))
        s["base_begintime"] = sbt
        s["base_record"] = tw
        Green.load!(s, env; showinfo=showinfo)
        detrendandtaper!(s["base_record"])
        detrendandtaper!(s["green_fun"])
        # taper!(s["green_fun"])
        for p in s["phases"], m in modules
            m.preprocess!(p, s, env)
        end
    end
    return nothing
end

"""
inverse!(env::Setting, modules::Vector{Module}, searchingMethod::Module)
-> (sdr, phaselist, misfit, misfitdetail)
"""
function inverse!(env::Setting, modules::Vector{Module}, searchingMethod::Module)
    phaselist = Tuple{Module,Setting}[]
    weightvec = Float64[]
    for s in env["stations"], m in modules
        for p in s["phases"]
            push!(phaselist, (m, p))
            push!(weightvec, m.weight(p, s, env))
        end
    end
    Lp = length(phaselist)
    sdr = Vector{Float64}[]
    misfitdetail = zeros(Float64, 0, Lp)
    misfit = Float64[]
    status = Setting()
    while searchingMethod.continueloop!(sdr, misfit, status, env)
        newsdr = searchingMethod.newparameters(sdr, misfit)
        newmomenttensor = dc2ts.(newsdr)
        newmisfitdetail = zeros(length(newsdr), Lp)
        newmisfit = zeros(Float64, length(newsdr))
        Threads.@threads for i = 1:(length(newsdr)*Lp)
            (p, q) = divrem(i - 1, Lp)
            p += 1
            q += 1
            func = phaselist[q][1].misfit
            @views newmisfitdetail[p, q] = func(phaselist[q][2], newmomenttensor[p])
        end
        append!(sdr, newsdr)
        misfitdetail = vcat(misfitdetail, newmisfitdetail)
        mul!(newmisfit, newmisfitdetail, weightvec)
        append!(misfit, newmisfit)
    end
    return (sdr, phaselist, misfit, misfitdetail)
end

function CAPmethod!(env::Setting, searchingMethod::Module)
    phaselist = Setting[]
    weightvec = Float64[]
    for s in env["stations"]
        for p in s["phases"]
            if !CAP.skip(p)
                push!(phaselist, p)
                push!(weightvec, CAP.weight(p, s, env))
            end
        end
    end
    Lp = length(phaselist)
    sdr = Vector{Float64}[]
    misfit = zeros(Float64, 0)
    sdr = searchingMethod.newparameters(sdr, misfit)
    Lm = length(sdr)
    momenttensor = dc2ts.(sdr)
    rec2 = zeros(Lp)
    syn2 = zeros(Lm, Lp)
    misfitdetail = zeros(Lm, Lp)
    misfit = zeros(Float64, Lm)
    Threads.@threads for q = 1:Lp
        @views rec2[q] = CAP.rec2(phaselist[q])
    end
    Threads.@threads for i = 1:(Lm*Lp)
        (p, q) = divrem(i - 1, Lp)
        p += 1
        q += 1
        @views syn2[p, q] = CAP.syn2(phaselist[q], momenttensor[p])
    end
    Lrec2 = sum(rec2)
    Lsyn2 = sum(syn2; dims = 2)
    Lsyn2 = reshape(Lsyn2, length(Lsyn2))
    M0 = sqrt.(Lrec2 ./ Lsyn2)

    Threads.@threads for i = 1:(Lm*Lp)
        (p, q) = divrem(i - 1, Lp)
        p += 1
        q += 1
        @views misfitdetail[p, q] = CAP.misfit(phaselist[q], momenttensor[p] .* M0[p])
    end
    mul!(misfit, misfitdetail, weightvec)
    return (sdr, phaselist, misfit, misfitdetail, M0)
end
