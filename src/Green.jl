"""
"""
module Green
using Printf, Dates, DelimitedFiles, Statistics, LinearAlgebra, SeismicRayTrace, Mmap, FFTW, DWN, SeisTools, TOML
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
    p = (t / tr - 8.0)
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
    # (r, baz, _) = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    r = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])*0.001
    baz = SeisTools.Geodesy.azimuth(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
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
# =   compressed 3D
# = =================

function _cglib_page_size(b::AbstractVector{<:Integer})
    p = zeros(Int, length(b))
    for i = eachindex(b)
        if i == 1
            p[i] = 1
        else
            p[i] = p[i-1] * b[i-1]
        end
    end
    return p
end

function _cglib_lin2cart(l::Int, b::AbstractVector{<:Integer})
    c = zeros(Int, length(b))
    p = _cglib_page_size(b)
    i = length(b)
    res = l - 1
    while i > 0
        (d, r) = divrem(res, p[i])
        c[i] = Int(d)
        res = r
        i -= 1
    end
    return c .+ 1
end

@inline function cglib_predict(x1::Real, x2::Real, x3::Real, c::Integer)
    if c == 0
        return Float64(x3)
    elseif c == 1
        return 2 * Float64(x3) - Float64(x2)
    elseif c == 2
        return Float64(x1) - 3 * Float64(x2) + 3 * Float64(x3)
    else
        return Float64(x2)
    end
end

@inline function cglib_decodenum(r::Unsigned, expmask::Unsigned,
    nsig::Integer, sigmask::Unsigned, expshift::Integer)
    a = r & sigmask
    b = (r >> nsig) & expmask
    e = Int(b) - expshift
    if a < 2^(nsig - 1)
        s = Float64(a) * 2.0^(e - nsig + 1)
    else
        s = (Float64(a) - 2.0^(nsig)) * 2.0^(e - nsig + 1)
    end
    return s
end

"""
cglib_readhead(io::IO) ->
(nx, ny, nz, x0, y0, z0, dx, dy, dz, nt, dt, stf,
grid,
bit_exp, expmask, expshift, bit_sig, sigmask, btype,
nleaf, leafv, nnode, leftnodes, rightnodes)
"""
function cglib_readhead(io::IO)
    flag = read(io, UInt8)
    partstart = zeros(UInt64, 4); read!(io, partstart)
    seek(io, 1 + 8 * 4 + 3 * 4)
    nt = read(io, Int32)
    dt = read(io, Float32)
    stf = zeros(Float32, nt); read!(io, stf)
    seek(io, partstart[1])
    ns = zeros(Int32, 3); read!(io, ns)
    xs = zeros(Float32, 6); read!(io, xs)
    (nx, ny, nz) = ns
    gridsize = Int.((nz, ny, nx))
    tracepospos = read(io, UInt64)
    seek(io, tracepospos)
    grid = zeros(UInt64, nz, ny, nx); read!(io, grid);

    seek(io, partstart[2])
    bit_exp = read(io, Int8)
    bit_sig = read(io, Int8)
    expshift = read(io, Int32)
    shiftval = read(io, Float64)
    nbit = bit_exp + bit_sig + 2

    typelist = (UInt8, UInt16, UInt32, UInt64)
    nbyte = ceil(Int, nextpow(2, nbit) / 8)
    btype = typelist[round(Int, log2(nbyte))+1]
    expmask = btype(2^bit_exp - 1)
    sigmask = btype(2^bit_sig - 1)

    nleaf = read(io, Int32)
    leafv = zeros(btype, nleaf); read!(io, leafv);
    nnode = read(io, Int32)
    leftnodes = zeros(Int32, nnode); read!(io, leftnodes);
    rightnodes = zeros(Int32, nnode); read!(io, rightnodes);

    return (nx, ny, nz, xs..., nt, dt, stf, grid,
        bit_exp, expmask, expshift, bit_sig, sigmask, btype,
        nleaf, leafv, nnode, leftnodes, rightnodes)
end

const cglib_bitorflag = (0b10000000,
             0b01000000,
             0b00100000,
             0b00010000,
             0b00001000,
             0b00000100,
             0b00000010,
             0b00000001);

function cglib_readtrace(io::IO, ix::Integer, iy::Integer, iz::Integer,
    nt::Integer, grid::Array{UInt64,3},
    nnode::Integer, leftnodes, rightnodes, leafv,
    bit_exp::Integer, expmask::Unsigned, expshift::Integer, bit_sig::Integer, sigmask::Unsigned
    )
    seek(io, grid[iz,iy,ix])
    tp = read(io, Float32)
    ts = read(io, Float32)
    nzero = read(io, Int32)
    decoded = zeros(Float32, nt, 6, 3);
    if nzero == nt
        return (tp, ts, decoded)
    end
    amp = read(io, Float32)
    cbyte = read(io, Int32)
    cbits = read(io, Int8)
    encoded = zeros(UInt8, cbyte); read!(io, encoded);

    tracedatasize = [nt-nzero, 6, 3]
    idata = 1
    ibyte = 0
    ibit = 8
    inode = nnode
    while true
        ibit += 1
        if ibit == 9
            ibit = 1
            ibyte += 1
        end
        inode = iszero(encoded[ibyte] & cglib_bitorflag[ibit]) ? leftnodes[inode] : rightnodes[inode]
        if iszero(leftnodes[inode])
            # println(idata)
            (it, im, ic) = _cglib_lin2cart(idata, tracedatasize)
            it += nzero
            pretype = (leafv[inode] >> (bit_exp+bit_sig)) & 0b11;
            if it == 1
                hp = 0.0
            elseif it == 2
                hp = cglib_predict(0.0, 0.0, decoded[it-1,im,ic], pretype)
            elseif it == 3
                hp = cglib_predict(0.0, decoded[it-2,im,ic], decoded[it-1,im,ic], pretype)
            else
                hp = cglib_predict(decoded[it-3,im,ic], decoded[it-2,im,ic], decoded[it-1,im,ic], pretype)
            end
            # decoded[it, im, ic] =  hp + dr * leafv[inode] - shiftval
            decoded[it,im,ic] = hp + cglib_decodenum(leafv[inode], expmask, bit_sig, sigmask, expshift)*Float64(amp);
            idata += 1
            inode = nnode
        end
        if ibyte == cbyte && ibit == cbits
            break
        end
    end
    return (tp, ts, decoded);
end

"""
cglib_readlocation(filename, x, y, z) -> (stf, dt, tp, ts, g)
"""
function cglib_readlocation(filename::AbstractString, x::Real, y::Real, z::Real)
    io = open(filename);
    (nx, ny, nz, x0, y0, z0, dx, dy, dz, nt, dt, stf, grid,
        bit_exp, expmask, expshift, bit_sig, sigmask, btype,
        nleaf, leafv, nnode, leftnodes, rightnodes) = cglib_readhead(io);
    if (x > (x0+(nx-1)*dx)) || (x < x0) ||
        (y > (y0+(ny-1)*dy)) || (y < y0) ||
        (z > (z0+(nz-1)*dz)) || (z < z0)
         error("Locaion out of range, require x($(x0),$(x0+(nx-1)*dx)), y($(y0),$(y0+(ny-1)*dy)), \
             z($(z0),$(z0+(nz-1)*dz)), current is x:$x, y:$y, z:$z")
     end
    xp = floor(Int, (x - x0) / dx) + 1;
    yp = floor(Int, (y - y0) / dy) + 1;
    zp = floor(Int, (z - z0) / dz) + 1;
    if xp == nx
        xp = nx-1;
    end
    if yp == ny
        yp = ny-1;
    end
    if zp == nz
        zp = nz - 1;
    end
    h = (x - x0) / dx - xp + 1.0
    k = (y - y0) / dy - yp + 1.0
    l = (z - z0) / dz - zp + 1.0
    parp = (nt, grid, nnode, leftnodes, rightnodes, leafv, bit_exp, expmask, expshift, bit_sig, sigmask)
    w = zeros(Float32, nt, 6, 3)
    gt1 = cglib_readtrace(io, xp,   yp,   zp,   parp...)
    gt2 = cglib_readtrace(io, xp+1, yp,   zp,   parp...)
    gt3 = cglib_readtrace(io, xp,   yp+1, zp,   parp...)
    gt4 = cglib_readtrace(io, xp+1, yp+1, zp,   parp...)
    gt5 = cglib_readtrace(io, xp,   yp,   zp+1, parp...)
    gt6 = cglib_readtrace(io, xp+1, yp,   zp+1, parp...)
    gt7 = cglib_readtrace(io, xp,   yp+1, zp+1, parp...)
    gt8 = cglib_readtrace(io, xp+1, yp+1, zp+1, parp...)
    w .+= gt1[3] .* ((1.0-h) * (1.0-k) * (1.0-l))
    w .+= gt2[3] .* (h       * (1.0-k) * (1.0-l))
    w .+= gt3[3] .* ((1.0-h) * k       * (1.0-l))
    w .+= gt4[3] .* (h       * k       * (1.0-l))
    w .+= gt5[3] .* ((1.0-h) * (1.0-k) * l)
    w .+= gt6[3] .* (h       * (1.0-k) * l)
    w .+= gt7[3] .* ((1.0-h) * k       * l)
    w .+= gt8[3] .* (h       * k       * l)
    tp = gt1[1] * ((1.0-h) * (1.0-k) * (1.0-l)) +
         gt2[1] * (h       * (1.0-k) * (1.0-l)) +
         gt3[1] * ((1.0-h) * k       * (1.0-l)) +
         gt4[1] * (h       * k       * (1.0-l)) +
         gt5[1] * ((1.0-h) * (1.0-k) * l) +
         gt6[1] * (h       * (1.0-k) * l) +
         gt7[1] * ((1.0-h) * k       * l) +
         gt8[1] * (h       * k       * l)
    ts = gt1[2] * ((1.0-h) * (1.0-k) * (1.0-l)) +
        gt2[2] * (h       * (1.0-k) * (1.0-l)) +
        gt3[2] * ((1.0-h) * k       * (1.0-l)) +
        gt4[2] * (h       * k       * (1.0-l)) +
        gt5[2] * ((1.0-h) * (1.0-k) * l) +
        gt6[2] * (h       * (1.0-k) * l) +
        gt7[2] * ((1.0-h) * k       * l) +
        gt8[2] * (h       * k       * l)
    close(io)
    return (stf, dt, tp, ts, w);
end


function load3dcompressedgreenlib(s, depth::Real, event, targetdir::AbstractString)
    # (r, baz, _) = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    r = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])*0.001
    baz = SeisTools.Geodesy.azimuth(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    x = r * cosd(baz)
    y = r * sind(baz)
    (_, dt, tp, ts, w) = cglib_readlocation(abspath(s["green_modelpath"]), x, y, depth)
    dt = Float64(dt)

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
                                "risetime"  => s["green_tsource"], # TODO bug here, this function should read risetime from glib
                                "tp"        => tp + 2.5 * s["green_tsource"],
                                "ts"        => ts + 2.5 * s["green_tsource"],
                                "type"      => s["green_modeltype"])
        printgreenfile(targetpath, meta, green[:, 6*(c-1).+(1:6)])
    end
end


# = =================
# =        2D
# = =================

function _glib_readlocation2d(filename::AbstractString, x::Real, z::Real)
    (n, xs, ys, zs, t, rt) = open(_glib_readhead, filename, "r")
    if (x > maximum(xs)) || (x < minimum(xs)) || (z > maximum(zs)) ||
       (z < minimum(zs))
        error("Locaion out of range in file $(filename),\nrequire x($(minimum(xs)),$(maximum(xs))), \
            z($(minimum(zs)),$(maximum(zs))), current is x:$x, z:$z")
    end
    ix = (x == xs[end]) ? length(xs) : max(2, findfirst(>(x), xs))
    iz = (z == zs[end]) ? length(zs) : max(2, findfirst(>(z), zs))
    w = zeros(Float32, Int(n[4]), 6, 3)
    io = open(filename, "r")
    H = Mmap.mmap(io, Array{Float32,6}, (Int(n[4]), 6, 3, Int(n[3]), Int(n[2]), Int(n[1])), Int((sum(n) + 5) * 4))
    h = (x - xs[ix-1]) / (xs[ix] - xs[ix-1])
    l = (z - zs[iz-1]) / (zs[iz] - zs[iz-1])
    for id = 1:3, ic = 1:6, it = 1:Int(n[4])
        w[it, ic, id] = H[it, ic, id, iz, 1, ix] * h * l +
                        H[it, ic, id, iz, 1, ix-1] * (1.0 - h) * l +
                        H[it, ic, id, iz-1, 1, ix] * h * (1.0 - l) +
                        H[it, ic, id, iz-1, 1, ix-1] * (1.0 - h) * (1.0 - l)
    end
    close(io)
    return (rt, t, w)
end

function load2dgreenlib(s, depth::Real, event, targetdir::AbstractString)
    # (r, baz, _) = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    r = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])*0.001
    baz = SeisTools.Geodesy.azimuth(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    raz = mod(baz+180.0, 360.0)
    (rt, t, w) = _glib_readlocation2d(abspath(s["green_modelpath"]), -r, depth)
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
    wz = -gt[:, 13:18]
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
# =   2D compressed
# = =================

function load2dcompressedgreenlib(s, depth::Real, event, targetdir::AbstractString)
    # (r, baz, _) = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    r = SeisTools.Geodesy.distance(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])*0.001
    baz = SeisTools.Geodesy.azimuth(s["meta_lat"], s["meta_lon"], event["latitude"], event["longitude"])
    raz = mod(baz+180.0, 360.0)
    (_, dt, tp, ts, w) = cglib_readlocation(abspath(s["green_modelpath"]), -r, 0.0, depth)
    setting = let
        (gfdir, gfilename) = splitdir(abspath(s["green_modelpath"]))
        TOML.parsefile(joinpath(gfdir, "setting.toml"))
    end
    rt = setting["risetime"]
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
    wz = -gt[:, 13:18]
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
    elseif uppercase(station["green_modeltype"]) == "3D_COMPRESSED"
        load3dcompressedgreenlib(station, env["algorithm"]["searchdepth"], env["event"], gfpath)
    elseif uppercase(station["green_modeltype"]) == "2D"
        load2dgreenlib(station, env["algorithm"]["searchdepth"], env["event"], gfpath)
    elseif uppercase(station["green_modeltype"]) == "2D_COMPRESSED"
        load2dcompressedgreenlib(station, env["algorithm"]["searchdepth"], env["event"], gfpath)
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
    @debug "load greenfun of station: $(station["network"]).$(station["station"])"
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
    @debug "green type: "*gmeta["type"]
    if gmeta["type"] == "DWN"
        (stf, _) = sourcetimefunction_v(npts, nfreq, gmeta["dt"] * npts, station["green_tsource"],
                                        -2 * station["green_tsource"], 1.0)
        g = zeros(size(tg))
        SeisTools.DataProcess.conv_f!(g, tg, stf)
    elseif uppercase(station["green_modeltype"]) == "3D"
        # (s1, _) = sourcetimefunction_v(npts, nfreq, gmeta["dt"]*npts, station["green_tsource"], -2*station["green_tsource"], 1.0)
        # S1 = fft(s1)
        # s2 = gauss.((0.0:npts-1).*gmeta["dt"], gmeta["risetime"])
        # S2 = fft(s2)
        # F = S1 .* conj.(S2) ./ max.(1e-5, abs2.(S2))
        # SeisTools.DataProcess.taper!(tg)
        # G = fft(tg)
        # for c in eachcol(G)
        #     c .*= F
        # end
        # ifft!(G)
        # g = real.(G)
        g = tg
    elseif uppercase(station["green_modeltype"]) == "3D_COMPRESSED"
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
    elseif uppercase(station["green_modeltype"]) == "2D"
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
    elseif uppercase(station["green_modeltype"]) == "2D_COMPRESSED"
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
    @debug "NPTS: $NPTS, Gnpts: $Gnpts, Nresample: $Nresample"
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
