"""
"""
module Green
using Printf, Dates, DelimitedFiles, Statistics, LinearAlgebra, SeismicRayTrace, Mmap, FFTW, DWN, SeisTools
import JuliaSourceMechanism: @must, @hadbetter, VelocityModel

properties = ["green_dt", "green_tsource", "green_model", "green_modeltype"]

_syn_lock = ReentrantLock()

# = ===============================
# =       Basic functions
# = ===============================
function _resample!(y::AbstractVector, x::AbstractVector)
    X = fft(x)
    n = length(y)
    Y = zeros(ComplexF64, n)
    m = floor(Int, min(n, length(x)) / 2)
    for i = 1:m
        Y[i] = X[i]
        Y[n-i+1] = X[end-i+1]
    end
    if mod(n, 2) > 0
        Y[m+1] = X[m+1]
    end
    ifft!(Y)
    for i = 1:n
        y[i] = real(Y[i])
    end
    return nothing
end

function _resample(x::AbstractVector, n::Int)
    @must mod(n, 2) == 0 "length must be n*2"
    y = zeros(n)
    _resample!(y, x)
    return y
end

function scalederf(t::Real, t0::Real)
    c = 6.0
    p = c * (t / t0 - 1.0)
    return (erf(p) + 1.0) / 2.0
end

function gauss(t::Real, t0::Real)
    tr = sqrt(-log(0.1)) * t0 / π
    p = (t / tr - 4.0)
    return exp(-p^2) / tr / sqrt(pi)
end

"""
"""
function intricker(t::Real, t0::Real)
    c = 6.0
    p = c * (t / t0 - 1.0)
    return -2 * c^2 * p * exp(-p^2) / (t0^2) / sqrt(pi)
end

function ricker(t::Real, t0::Real)
    c = 6.0
    p = c * (t / t0 - 1.0)
    return -2 * c^3 * (1.0 - 2.0 * p^2) * exp(-p^2) / (t0^3) / sqrt(pi)
end

# = =================
# =        DWN
# = =================

"""
_calgreenfun_dwn_station(station, model::Matrix, depth::Real, event, strs) -> nothing

input:

  - station: station Dict
  - model: Matrix{Float} with format `depth,Vp,Vs,density,Qp,Qs`
  - depth: depth of source
  - event: Dict
  - strs: NamedTuple, paths that are used, including tmp(temporary dir), model(model name), green(dir to store result)
"""
function _calgreenfun_dwn_station(s, model::Matrix, depth::Real, event, strs)
    # initial
    zr = 0.0
    el = round(s["meta_el"] / 1000.0, digits=3)
    zh = round(depth + el, digits=3)
    rd = s["base_distance"]
    az = s["base_azimuth"]
    dep = round.(model[:, 1] .+ el, digits=3)
    # process model
    l = findall(dep .> 0.0)
    if l[1] > 1
        l = [l[1] - 1; l]
    end
    m = model[l, :]
    m[:, 1] = m[:, 1] .+ el
    m[1, 1] = 0.0
    m[:, 1] = [diff(m[:, 1]); 0.0]
    if m[1, 1] < 0.01
        m = m[2:end, :]
    end
    m = round.(m, digits = 6)
    npts = round(Int, s["green_tl"] / s["green_dt"])
    npts += mod(npts, 2)
    if "green_xl" in keys(s)
        xl = s["green_xl"]
    else
        xl = ceil((s["green_dt"] * npts * maximum(m[:, 2]) + s["base_distance"]) / 10.0) * 10.0
        s["green_xl"] = xl
    end
    spec = DWN.dwn(m, zh, 1.0, [(rd, az)], zr, npts, s["green_dt"], xl, s["green_dt"], s["green_m"])
    w = DWN.freqspec2timeseries(spec, npts)
    green = [w[1, 2, :]... w[1, 1, :]... -w[1, 3, :]...]
    targetdir = strs.green
    mkpath(targetdir)
    cmps = ["E", "N", "Z"]
    (tp, ts) = let
        m0 = deepcopy(m)
        m0[2:end, 1] .= cumsum(m0[:, 1])[1:end-1]
        m0[1, 1] = 0.0
        m0 = round.(m0, digits=6)
        tp = raytrace_fastest(0.0, zh, rd, m0[:, 1], m0[:, 2])
        ts = raytrace_fastest(0.0, zh, rd, m0[:, 1], m0[:, 3])
        (tp.phase.t, ts.phase.t)
    end
    for c = 1:3
        targetpath = joinpath(targetdir, @sprintf("%s.%s.%s.gf", s["network"], s["station"], cmps[c]))
        meta = Dict{String,Any}("modelname" => strs.model,
                                "network"   => s["network"],
                                "station"   => s["station"],
                                "component" => s["component"],
                                "depth"     => zh,
                                "distance"  => rd,
                                "azimuth"   => az,
                                "dt"        => s["green_dt"],
                                "model"     => m,
                                "tp"        => tp,
                                "ts"        => ts,
                                "type"      => "DWN")
        printgreenfile(targetpath, meta, green[:, 6*(c-1).+(1:6)])
    end
    return nothing
end

"""
calgreenfun_dwn(station, model::Matrix, depth::Real, event, strs)

a front end of `_calgreenfun_dwn_station`, if the Greenfun is already exist, skip the station.
! will be edit to support more feature
"""
function calgreenfun_dwn(station, model::Matrix, depth::Real, event, strs)
    targetfilenames = map(c -> abspath(strs.green, @sprintf("%s.%s.%s.gf", station["network"], station["station"], c)),
                          ["E", "N", "Z"])
    if !all(isfile.(targetfilenames))
        _calgreenfun_dwn_station(station, model, depth, event, strs)
    end
    return nothing
end

# = =================
# =      SEM/FD
# = =================
function _glib_readhead(io::IO)
    rt = read(io, Float32)
    n = zeros(Int32, 4)
    read!(io, n)
    x = zeros(Float32, n[1])
    y = zeros(Float32, n[2])
    z = zeros(Float32, n[3])
    t = zeros(Float32, n[4])
    read!(io, x)
    read!(io, y)
    read!(io, z)
    read!(io, t)
    return (n = n, x = x, y = y, z = z, t = t, risetime = rt)
end

function _glib_readall(io::IO)
    (_, x, y, z, t, rt) = _glib_readhead(io)
    H = zeros(Float32, length(t), 6, 3, length(z), length(y), length(x))
    read!(io, H)
    return (x, y, z, t, H, rt)
end

function _glib_readall(s::AbstractString)
    return open(_glib_readall, s, "r")
end

function _glib_readlocation(filename::AbstractString, x::Real, y::Real, z::Real)
    (n, xs, ys, zs, t, rt) = open(_glib_readhead, filename, "r")
    if (x > maximum(xs)) || (x < minimum(xs)) || (y > maximum(ys)) || (y < minimum(ys)) || (z > maximum(zs)) ||
       (z < minimum(zs))
        error("Locaion out of range in file $(filename),\nrequire x($(minimum(xs)),$(maximum(xs))), y($(minimum(ys)),$(maximum(ys))), \
            z($(minimum(zs)),$(maximum(zs))), current is x:$x, y:$y, z:$z")
    end
    ix = (x == xs[end]) ? length(xs) : max(2, findfirst(>(x), xs))
    iy = (y == ys[end]) ? length(ys) : max(2, findfirst(>(y), ys))
    iz = (z == zs[end]) ? length(zs) : max(2, findfirst(>(z), zs))
    w = zeros(Float32, Int(n[4]), 6, 3)
    io = open(filename, "r")
    H = Mmap.mmap(io, Array{Float32,6}, (Int(n[4]), 6, 3, Int(n[3]), Int(n[2]), Int(n[1])), Int((sum(n) + 5) * 4))
    h = (x - xs[ix-1]) / (xs[ix] - xs[ix-1])
    k = (y - ys[iy-1]) / (ys[iy] - ys[iy-1])
    l = (z - zs[iz-1]) / (zs[iz] - zs[iz-1])
    for id = 1:3, ic = 1:6, it = 1:Int(n[4])
        w[it, ic, id] = H[it, ic, id, iz, iy, ix] * h * k * l +
                        H[it, ic, id, iz, iy, ix-1] * (1.0 - h) * k * l +
                        H[it, ic, id, iz, iy-1, ix] * h * (1.0 - k) * l +
                        H[it, ic, id, iz, iy-1, ix-1] * (1.0 - h) * (1.0 - k) * l +
                        H[it, ic, id, iz-1, iy, ix] * h * k * (1.0 - l) +
                        H[it, ic, id, iz-1, iy, ix-1] * (1.0 - h) * k * (1.0 - l) +
                        H[it, ic, id, iz-1, iy-1, ix] * h * (1.0 - k) * (1.0 - l) +
                        H[it, ic, id, iz-1, iy-1, ix-1] * (1.0 - h) * (1.0 - k) * (1.0 - l)
    end
    close(io)
    return (rt, t, w)
end

"""
ttlib_readhead(io::IO) -> (nx, ny, nz, dx, dy, dz, ox, oy, oz)
"""
function ttlib_readhead(io::IO)
    (nx, ny, nz) = _read_vector(io, Int32, 3)
    (dx, dy, dz, ox, oy, oz) = _read_vector(io, Float32, 6)
    return (nx, ny, nz, dx, dy, dz, ox, oy, oz)
end

function _read_vector(io::IO, T::Type, n::Integer)
    t = zeros(T, n)
    read!(io, t)
    return t
end

"""
ttlib_readall(io::IO) -> (xs, ys, zs, tt)
"""
function ttlib_readall(io::IO)
    (nx, ny, nz, dx, dy, dz, ox, oy, oz) = ttlib_readhead(io)
    tt = zeros(Float32, 2, nz, ny, nx)
    read!(io, tt)
    return (dx, dy, dz, ox, oy, oz, tt)
end

"""
ttlib_readlocation(io::IO, x::Real, y::Real, z::Real) -> (tp, ts)
"""
function ttlib_readlocation(io::IO, x::Real, y::Real, z::Real)
    (nx, ny, nz, dx, dy, dz, ox, oy, oz) = ttlib_readhead(io)
    ix = floor(Int, (x - ox) / dx) + 1
    iy = floor(Int, (y - oy) / dy) + 1
    iz = floor(Int, (z - oz) / dz) + 1
    if (ix < 1) || (iy < 1) || (iz < 1) || (ix > nx) || (iy > ny) || (iz > nz)
        error("Location out of range")
    end
    if ix == nx
        ix -= 1
    end
    if iy == ny
        iy -= 1
    end
    if iz == nz
        iz -= 1
    end
    h = (x - ox) / dx - ix + 1.0
    k = (y - oy) / dy - iy + 1.0
    l = (z - oz) / dz - iz + 1.0
    seek(io, ((ix - 1) * 2 * nz * ny + (iy - 1) * 2 * nz + (iz - 1) * 2 + 1 - 1 + 9) * 4)
    (p000, s000) = _read_vector(io, Float32, 2)
    seek(io, (ix * 2 * nz * ny + (iy - 1) * 2 * nz + (iz - 1) * 2 + 1 - 1 + 9) * 4)
    (p100, s100) = _read_vector(io, Float32, 2)
    seek(io, ((ix - 1) * 2 * nz * ny + iy * 2 * nz + (iz - 1) * 2 + 1 - 1 + 9) * 4)
    (p010, s010) = _read_vector(io, Float32, 2)
    seek(io, (ix * 2 * nz * ny + iy * 2 * nz + (iz - 1) * 2 + 1 - 1 + 9) * 4)
    (p110, s110) = _read_vector(io, Float32, 2)
    seek(io, ((ix - 1) * 2 * nz * ny + (iy - 1) * 2 * nz + iz * 2 + 1 - 1 + 9) * 4)
    (p001, s001) = _read_vector(io, Float32, 2)
    seek(io, (ix * 2 * nz * ny + (iy - 1) * 2 * nz + iz * 2 + 1 - 1 + 9) * 4)
    (p101, s101) = _read_vector(io, Float32, 2)
    seek(io, ((ix - 1) * 2 * nz * ny + iy * 2 * nz + iz * 2 + 1 - 1 + 9) * 4)
    (p011, s011) = _read_vector(io, Float32, 2)
    seek(io, (ix * 2 * nz * ny + iy * 2 * nz + iz * 2 + 1 - 1 + 9) * 4)
    (p111, s111) = _read_vector(io, Float32, 2)
    tp = p000 * (1.0 - h) * (1.0 - k) * (1.0 - l) +
         p100 * h * (1.0 - k) * (1.0 - l) +
         p010 * (1.0 - h) * k * (1.0 - l) +
         p110 * h * k * (1.0 - l) +
         p001 * (1.0 - h) * (1.0 - k) * l +
         p101 * h * (1.0 - k) * l +
         p011 * (1.0 - h) * k * l +
         p111 * h * k * l
    ts = s000 * (1.0 - h) * (1.0 - k) * (1.0 - l) +
         s100 * h * (1.0 - k) * (1.0 - l) +
         s010 * (1.0 - h) * k * (1.0 - l) +
         s110 * h * k * (1.0 - l) +
         s001 * (1.0 - h) * (1.0 - k) * l +
         s101 * h * (1.0 - k) * l +
         s011 * (1.0 - h) * k * l +
         s111 * h * k * l
    return (tp, ts)
end

function ttlib_readlocation(path::AbstractString, x::Real, y::Real, z::Real)
    io = open(path)
    (tp, ts) = ttlib_readlocation(io, x, y, z)
    close(io)
    return (tp, ts)
end

function load3dgreenlib(s, depth::Real, event, targetdir::AbstractString)
    (r, baz, _) = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    x = r * cosd(baz)
    y = r * sind(baz)
    (rt, t, w) = _glib_readlocation(abspath(s["green_modelpath"]), x, y, depth)
    dt = t[2] - t[1] |> Float64

    # (pidx, sidx) = autopick(w[:, :, 2], w[:, :, 1], s["base_azimuth"]; LW = round(Int, 10 * rt / dt),
    #                         SW = round(Int, 3 * rt / dt))
    # tp = pidx * dt
    # ts = sidx * dt
    (tp, ts) = ttlib_readlocation(abspath(s["green_ttlibpath"]), x, y, depth)
    mkpath(targetdir)
    cmps = ["E", "N", "Z"]
    green = [w[:, :, 2] w[:, :, 1] -w[:, :, 3]]
    for c = 1:3
        targetpath = joinpath(targetdir, @sprintf("%s.%s.%s.gf", s["network"], s["station"], cmps[c]))
        meta = Dict{String,Any}("modelname" => s["green_model"],
                                "network"   => s["network"],
                                "station"   => s["station"],
                                "component" => s["component"],
                                "depth"     => depth,
                                "distance"  => r,
                                "bazimuth"  => baz,
                                "dt"        => dt,
                                "risetime"  => rt,
                                "tp"        => tp + 0.5 * rt,
                                "ts"        => ts + 0.5 * rt,
                                "type"      => s["green_modeltype"])
        printgreenfile(targetpath, meta, green[:, 6*(c-1).+(1:6)])
    end
end

# = =================
# =        2D
# = =================

function load2dgreenlib(s, depth::Real, event, targetdir::AbstractString)
    (r, baz, _) = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    raz = mod(baz+180.0, 360.0)
    (rt, t, w) = _glib_readlocation(abspath(s["green_modelpath"]), -r, 0.0, depth)
    dt = t[2] - t[1] |> Float64
    (tp, ts) = ttlib_readlocation(abspath(s["green_ttlibpath"]), -r, 0.0, depth)
    mkpath(targetdir)
    rrmat=[1 4 5; 4 2 6; 5 6 3]
    T = [cosd(raz) sind(raz) 0.0;
        -sind(raz) cosd(raz) 0.0;
        0.0 0.0 1.0]
    Tt = zeros(18, 18)
    for i = 1:3, j = 1:3, n = 1:3, p = 1:3, q = 1:3, m = 1:3
        r = rrmat[i, j] + (n-1)*6
        c = rrmat[p, q] + (m-1)*6
        f = rrmat[p, q] > 3 ? 0.5 : 1.0
        Tt[c, r] += T[m, n]*T[p, i]*T[q, j]*f
    end
    wr = w[:, :, 1]
    wt = w[:, :, 2]
    wd = w[:, :, 3]
    gt = [wr wt wd] * Tt
    we = gt[:, 7:12]
    wn = gt[:, 1:6]
    wz = -wd[:, 13:18]
    cmps = ["E", "N", "Z"]
    green = [we wn wz]
    for c = 1:3
        targetpath = joinpath(targetdir, @sprintf("%s.%s.%s.gf", s["network"], s["station"], cmps[c]))
        meta = Dict{String,Any}("modelname" => s["green_model"],
                                "network"   => s["network"],
                                "station"   => s["station"],
                                "component" => s["component"],
                                "depth"     => depth,
                                "distance"  => r,
                                "bazimuth"  => baz,
                                "dt"        => dt,
                                "risetime"  => rt,
                                "tp"        => tp + 0.5 * rt,
                                "ts"        => ts + 0.5 * rt,
                                "type"      => s["green_modeltype"])
        printgreenfile(targetpath, meta, green[:, 6*(c-1).+(1:6)])
    end
end

# = =================
# =        IO
# = =================
function mat2line(x::AbstractMatrix, intp::AbstractString)
    return map(1:size(x, 1)) do i
        join(x[i, :], intp)
    end
end

function scangreenfile(path::AbstractString)
    meta = Dict{String,Any}()
    for l in filter(v -> contains(v, "#"), readlines(path))
        sl = split(l, " "; keepempty = false)
        k = String(sl[2])
        t = eval(Meta.parse(sl[3]))
        if k == "layer"
            k = "model"
            if "model" in keys(meta)
                v = vcat(meta["model"], permutedims(parse.(Float64, sl[3:end])))
            else
                v = permutedims(parse.(Float64, sl[3:end]))
            end
        elseif k == "model"
            v = zeros(0, 6)
        elseif t <: Real
            v = parse(t, sl[4])
        else
            v = String(join(sl[4:end], ' '))
        end
        meta[k] = v
    end
    greendata = readdlm(path, ','; comments = true)
    return (meta, greendata)
end

function printgreenfile(path::AbstractString, meta::Dict{String,Any}, gf::AbstractMatrix)
    if isfile(path)
        @error "file already exist"
    end
    open(path, "w") do io
        for k in sort(collect(keys(meta)))
            if typeof(meta[k]) <: AbstractMatrix
                s = "# " * string(k) * " Matrix\n" * join("# layer " .* mat2line(meta[k], " "), "\n")
            else
                s = "# " * string(k) * " " * string(typeof(meta[k])) * " " * string(meta[k])
            end
            println(io, s)
        end
        for l in mat2line(gf, ",")
            println(io, l)
        end
    end
    return nothing
end

function _conv_timedomain(u::Matrix, v::Vector)
    c = zeros(size(u))
    Lv = length(v)
    # for i = 1:size(c, 1), j = 1:size(c, 2)
    #     for k = 1:size(c, 1)
    #         if i - k < 0 || i - k + 1 > Lv
    #             continue
    #         end
    #         c[i, j] += u[k, j] * v[i-k+1]
    #     end
    # end
    for j in axes(c, 2), i in axes(c, 1), k in axes(c, 1)
        c[i, j] += u[k, j] * v[mod(i - k, Lv)+1]
    end
    return c
end

function _conv_freqdomain(u::Matrix, V::Vector)
    w = zeros(size(u))
    for i = 1:size(u, 2)
        U = fft(u[:, i])
        w[:, i] .= ifft(U .* V) .|> real
    end
    return w
end

function greenfilename(station::Dict, env::Dict)
    gfpath = normpath(env["dataroot"], "greenfun",
                      @sprintf("%s-%.4f",
                               station["green_model"], env["algorithm"]["searchdepth"]))
    gfname = station["network"] * "." * station["station"] * "." * station["component"] * ".gf"
    return (gfpath, gfname)
end

function calculategreenfun(station::Dict, env::Dict)
    (gfpath, _) = greenfilename(station, env)
    if uppercase(station["green_modeltype"]) == "DWN"
        modelpath = normpath(env["dataroot"], "model", station["green_model"] * ".model")
        if !isfile(modelpath)
            if station["green_model"] == "crust1.0"
                lock(_syn_lock)
                try
                    (modeldir, _) = splitdir(modelpath)
                    if !isdir(modeldir)
                        mkpath(modeldir)
                    end
                    _m = VelocityModel.readmodel_crust10(env["event"]["latitude"], env["event"]["longitude"])
                    _t = zeros(size(_m, 1), 6)
                    _t[:, 1:4] .= _m
                    _t[:, 5] .= 200.0
                    _t[:, 6] .= 100.0
                    writedlm(modelpath, _t, ',')
                catch _err
                    throw(_err)
                finally
                    unlock(_syn_lock)
                end
            else
                error("Model file not exist", station["network"] * "." * station["station"])
            end
        end
        model = readdlm(modelpath, ','; comments = true)
        calgreenfun_dwn(station, model, env["algorithm"]["searchdepth"], env["event"],
                        (model = station["green_model"], green = gfpath))
    elseif uppercase(station["green_modeltype"]) == "3D"
        load3dgreenlib(station, env["algorithm"]["searchdepth"], env["event"], gfpath)
    else
        error("Unknown Green lib type, station: ", station["network"] * "." * station["station"])
    end
    return nothing
end

"""
load(station::Dict, env::Dict) -> g

load Green function
"""
function load!(station::Dict, env::Dict; showinfo::Bool=false)
    (gfpath, gfname) = greenfilename(station, env)
    gfilename = joinpath(gfpath, gfname)
    if !isfile(gfilename)
        if showinfo
            @info "Green function of $(station["network"]).$(station["station"]).$(station["component"]) not exist. Calculat now"
        end
        calculategreenfun(station, env)
    end
    (gmeta, tg) = scangreenfile(gfilename)
    npts = size(tg, 1)
    nfreq = round(Int, npts / 2)
    if gmeta["type"] == "DWN"
        (stf, _) = sourcetimefunction_v(npts, nfreq, gmeta["dt"] * npts, station["green_tsource"],
                                        -2 * station["green_tsource"], 1.0)
        g = zeros(size(tg))
        SeisTools.DataProcess.conv_f!(g, tg, stf)
    elseif uppercase(station["green_modeltype"]) == "3D"
        (s1, _) = sourcetimefunction_v(npts, nfreq, gmeta["dt"]*npts, station["green_tsource"], -2*station["green_tsource"], 1.0)
        S1 = fft(s1)
        s2 = gauss.((0.0:npts-1).*gmeta["dt"], gmeta["risetime"])
        S2 = fft(s2)
        F = S1 .* conj.(S2) ./ max.(1e-5, abs2.(S2))
        SeisTools.DataProcess.taper!(tg)
        G = fft(tg)
        for c in eachcol(G)
            c .*= F
        end
        ifft!(G)
        g = real.(G)
    end
    shift = round(Int,
                  Millisecond(env["event"]["origintime"] - station["base_trim"][1]) /
                  Millisecond(round(Int, gmeta["dt"] * 1000)))
    NPTS = round(Int,
                 Millisecond(station["base_trim"][2] - station["base_trim"][1]) /
                 Millisecond(round(Int, gmeta["dt"] * 1000)))
    Gnpts = min(round(Int,
                      Millisecond(station["base_trim"][2] - env["event"]["origintime"]) /
                      Millisecond(round(Int, gmeta["dt"] * 1000))), npts)
    Nresample = round(Int, Gnpts * gmeta["dt"] / station["green_dt"])
    if shift < 0
        @error("Station: " * station["network"] * "." * station["station"] * " shift is less than 0: " * string(shift))
    end
    if Nresample + shift <= 0
        error("Station: " * station["network"] * "." * station["station"] * " length of Green function too short")
    end
    wr = zeros(NPTS, 6)
    for i = 1:6
        # _resample!(@view(wr[shift+1:shift+Nresample, i]), g[:, i])
        SeisTools.DataProcess.resample!(@view(wr[max(1, shift + 1):shift+Nresample, i]), g[1:Gnpts, i])
    end
    station["green_fun"] = wr
    station["green_dt"] = gmeta["dt"]
    return nothing
end

function calc(station::Dict, env::Dict; showinfo::Bool = false)
    (gfpath, gfname) = greenfilename(station, env)
    gfilename = joinpath(gfpath, gfname)
    if !isfile(gfilename)
        if showinfo
            @info "Green function of $(station["network"]).$(station["station"]).$(station["component"]) not exist. Calculat now"
        end
        calculategreenfun(station, env)
    end
    return nothing
end

end # module Green
