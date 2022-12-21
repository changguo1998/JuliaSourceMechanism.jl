#=
"""
distance(lat1, lon1, lat2, lon2) -> (dist, az, gcarc)

compute distance and azimuth between two points on the Earth according to the reference Sphere.
distance is in km, az in degree, centered at point 1, and garc in radius degree
"""
function distance(lat1, lon1, lat2, lon2)
    R = 6371.0
    θ1 = deg2rad(90.0 - lat1)
    θ2 = deg2rad(90.0 - lat2)
    φ1 = deg2rad(lon1)
    φ2 = deg2rad(lon2)
    n1 = [sin(θ1) * cos(φ1), sin(θ1) * sin(φ1), cos(θ1)]
    n2 = [sin(θ2) * cos(φ2), sin(θ2) * sin(φ2), cos(θ2)]
    gcarc = acos(n1' * n2)
    dist = R * gcarc
    t12 = normalize(n2 .- (n1' * n2) .* n1)
    tnorth = [sin(θ1 - pi / 2) * cos(φ1), sin(θ1 - pi / 2) * sin(φ1), cos(θ1 - pi / 2)]
    teast = [cos(φ1 + pi / 2), sin(φ1 + pi / 2), 0.0]
    az = atand(t12' * teast, t12' * tnorth) |> x -> mod(x, 360.0)
    return (dist, az, gcarc)
end

function taper!(x::VecOrMat; ratio::Real = 0.05)
    # taper!(t->t, x, ratio=ratio)
    taper!(t -> 1.0 - cos(t * pi / 2), x; ratio = ratio)
    return nothing
end

function taper!(f::Function, x::VecOrMat; ratio::Real = 0.05)
    N = size(x, 1)
    M = max(round(Int, N * ratio), 2)
    for i = 1:M, j = axes(x, 2)
        x[i, j] *= f((i - 1) / (M - 1))
        x[N-i+1, j] *= f((i - 1) / (M - 1))
    end
    return nothing
end

"""
detrendandtaper!(x; ratio::Real=0.05)

remove least square linear trend and then taper using cosine window
"""
function detrendandtaper!(x::VecOrMat; ratio::Real = 0.05)
    N = size(x, 1)
    ibar = (1 + N) / 2.0
    xbar = mean(x; dims = 1)[:]
    si = N * (N * N - 1) / 12
    crossm = zeros(size(x, 2))
    for i = axes(x, 1), j = axes(x, 2)
        crossm[j] += i * x[i, j]
    end
    k = @. (crossm - ibar * xbar * N) / si
    b = @. xbar - k * ibar
    for i = 1:N, j = axes(x, 2)
        x[i, j] -= k[j] * i + b[j]
    end
    taper!(x; ratio = ratio)
    return nothing
end
=#

function pick_windowratio(x::AbstractVector{<:Real}, wl::Integer)
    L = length(x)
    r = zeros(L - 2*wl)
    for i = eachindex(r)
        r[i] = std(x[i+wl:i+2*wl])/std(x[i:i+wl])
    end
    (_, j) = findmax(r)
    return (j+wl, r)
end

function _freedom(x::AbstractVecOrMat)
    if size(x, 2) == 1
        return 1.0
    else
        F = svd(x)
        return sum(F.S)/maximum(F.S)
    end
end

function pick_freedom(x::AbstractVecOrMat{<:Real}, wl::Integer)
    L = size(x, 1)
    r = zeros(L-wl)
    for i = eachindex(r)
        r[i] = _freedom(x[i:i+wl, :])
    end
    (_, j) = findmin(r)
    return (j+wl, r)
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

#=
function sacDateTime(h)
    return DateTime(h["nzyear"], 1, 1, h["nzhour"], h["nzmin"], h["nzsec"], h["nzmsec"]) + Day(h["nzjday"] - 1)
end

function trim(w::VecOrMat, wbt::Real, bt::Real, et::Real, dt::Real; fillval::Union{Nothing,Real} = nothing)
    npts = round(Int, (et - bt) / dt)
    bi = round(Int, (bt - wbt) / dt) + 1
    if isnothing(fillval)
        ei = bi + npts - 1
        bi = max(1, bi)
        ei = min(ei, length(w))
        return (ndims(w) == 1) ? w[bi:ei] : w[bi:ei, :]
    else
        if ndims(w) == 1
            nw = fill(fillval, npts)
        else
            nw = fill(fillval, npts, size(w, 2))
        end
        for i = 1:npts, j = axes(w, 2)
            if (bi + i - 1 > 0) && (bi + i - 1 < size(w, 1))
                nw[i, j] = w[bi+i-1, j]
            end
        end
        return nw
    end
end

function trim(w::VecOrMat, wbt::DateTime, bt::DateTime, et::DateTime, dt::Real; fillval::Union{Nothing,Real} = nothing)
    return trim(w, 0.0, round(bt - wbt, Millisecond).value * 1e-3, round(et - wbt, Millisecond).value * 1e-3, dt;
                fillval = fillval)
end

function trim(sac::Seis.SACFrame, bt::DateTime, et::DateTime; fillval::Union{Nothing,Real} = nothing)
    nh = deepcopy(sac.head)
    npts = round(Int, round(et - bt, Millisecond) / Millisecond(round(Int, sac.head["delta"] * 1000)))
    shift = round(Int,
                  round(bt - sacDateTime(sac.head), Millisecond) / Millisecond(round(Int, sac.head["delta"] * 1000))) +
            1
    if isnothing(fillval)
        esh = shift - 1 + npts
        shift = max(1, shift)
        esh = min(length(sac.data[1]), esh)
        w = deepcopy(sac.data[1][shift:esh])
    else
        w = fill(fillval, npts)
        for i = 1:npts
            if (shift - 1 + i > 1) && (shift - 1 + i < length(sac.data[1]))
                w[i] = sac.data[1][shift-1+i]
            end
        end
    end
    nh["nzyear"] = year(bt)
    nh["nzjday"] = dayofyear(bt)
    nh["nzhour"] = hour(bt)
    nh["nzmin"] = minute(bt)
    nh["nzsec"] = second(bt)
    nh["nzmsec"] = millisecond(bt)
    nh["npts"] = length(w)
    nh["depmen"] = mean(w)
    nh["depmax"] = maximum(w)
    nh["depmin"] = minimum(w)
    nh["b"] = 0.0
    return Seis.SACFrame(nh, [w])
end

function trim(s::Setting, bt::DateTime, et::DateTime; fillval::Union{Nothing,Real} = nothing)
    return trim(s["base_record"], s["base_begintime"], bt, et, s["meta_dt"]; fillval = fillval)
end

=#

function calcgreen!(env::Setting)
    taglist = String[]
    idxlist = Int[]
    for i in eachindex(env["stations"])
        s = env["stations"][i]
        tag = String(s["network"]*"."*s["station"])
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
        Green.calc(s, env)
    end
    return nothing
end

function loaddata!(env::Setting)
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
        Green.load!(s, env)
        detrendandtaper!(s["base_record"])
        detrendandtaper!(s["green_fun"])
    end
    return nothing
end

"""
preprocess!(env::Setting, modules::Vector{Module})
"""
function preprocess!(env::Setting, modules::Vector{Module})
    for s in env["stations"]
        t = SeisTools.SAC.read(normpath(env["dataroot"], "sac", s["meta_file"]))
        trim_bt = s["base_trim"][1]
        trim_et = s["base_trim"][2]
        (sbt, tw, _) = SeisTools.DataProcess.cut(t.data, SeisTools.SAC.DateTime(t.hdr), trim_bt, trim_et, 
            Millisecond(round(Int, t.hdr["delta"]*1000)))
        s["base_begintime"] = sbt
        s["base_record"] = tw
        Green.load!(s, env)
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
