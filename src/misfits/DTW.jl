"""
JINV_dtw

misfit function using DTW algorithm
"""
module DTW

using Dates, DSP, Statistics, LinearAlgebra
import JuliaSourceMechanism: trim, Setting

tags = ("dtw", "DTW")
properties = ["dtw_dt", "dtw_maxlag", "dtw_klim", "dtw_trim", "dtw_order", "dtw_band"]

function weight(p::Setting, s::Setting, e::Setting)
    if "dtw_weight" ∈ keys(p)
        return p["dtw_weight"]
    else
        idx = findfirst(x -> x ∈ tags, e["algorithm"]["misfit"])
        if isnothing(idx)
            return 1.0
        else
            return e["algorithm"]["weight"][idx]
        end
    end
end

function preprocess!(phase::Setting, station::Setting, env::Setting)
    w = deepcopy(station["base_record"])
    w_resample = resample(w, station["meta_dt"] / phase["dtw_dt"])
    fltr = digitalfilter(Bandpass(phase["dtw_band"][1], phase["dtw_band"][2]; fs = 1 / phase["dtw_dt"]),
                         Butterworth(phase["dtw_order"]))
    w_filt = filtfilt(fltr, w_resample)
    w_trim = trim(w_filt, station["base_begintime"],
                  phase["at"] + Millisecond(round(Int, phase["dtw_trim"][1] * 1e3)),
                  phase["at"] + Millisecond(round(Int, phase["dtw_trim"][2] * 1e3)), phase["dtw_dt"])
    nm = norm(w_trim)
    w_trim ./= nm
    g = deepcopy(station["green_fun"])
    g_resample = resample(g, station["green_dt"] / phase["dtw_dt"]; dims = 1)
    g_filt = filtfilt(fltr, g_resample)
    g_trim = trim(g_filt, env["event"]["origintime"],
                  env["event"]["origintime"] + Millisecond(round(Int, (phase["tt"] + phase["dtw_trim"][1]) * 1e3)),
                  env["event"]["origintime"] + Millisecond(round(Int, (phase["tt"] + phase["dtw_trim"][2]) * 1e3)),
                  phase["dtw_dt"]; fillval = 0.0)
    phase["dtw_record"] = w_trim
    phase["dtw_greenfun"] = g_trim
    phase["dtw_Imaxlag"] = round(Int, phase["dtw_maxlag"] / phase["dtw_dt"])
    phase["dtw_Iklim"] = round(Int, phase["dtw_klim"] / phase["dtw_dt"])
    return nothing
end

function errormap(x::Vector, g::Matrix, m::Vector, maxlag::Int)
    N = length(x)
    e = zeros(N, 2 * maxlag + 1)
    w = normalize(g * m)
    for i = 1:N, j = 1:2*maxlag+1
        iw = min(max(1, i + j - maxlag - 1), N)
        e[i, j] = (x[i] - w[iw])^2
    end
    return e
end

function cumulate(e::Matrix, klim::Int)
    maxlag = round(Int, (size(e, 2) - 1) / 2)
    c = zeros(typeof(e[1, 1]), size(e))
    for i = 1:size(e, 1)
        imax = max(1, min(size(e, 1), i - 1))
        imin = max(1, min(size(e, 1), i - klim))
        for l = -maxlag:maxlag
            j = min(max(1, l + maxlag + 1), 2 * maxlag + 1)
            lUp = min(2 * maxlag + 1, j + 1)
            lDown = max(1, j - 1)
            cUp = 0.0
            cDown = 0.0
            for erri = imin+1:imax
                cUp += e[erri, lUp]
                cDown += e[erri, lDown]
            end
            cUp += c[imin, lUp]
            cDown += c[imin, lDown]
            cMiddle = c[imax, j]
            c[i, j] = e[i, j] + min(cUp, cDown, cMiddle)
        end
    end
    return c
end

function trace(e::Matrix{T}, c::Matrix{T}, klim::Int) where {T<:Real}
    (N, L) = size(c)
    maxlag = round(Int, (L - 1) / 2)
    p1 = zeros(Int, N)
    p2 = zeros(Int, N)
    k = N
    t1 = findmin(c[N, :])
    p1[k] = N
    p2[k] = t1[2]
    # p2[k] = maxlag + 1

    k -= 1
    while k > 0
        i = p1[k+1]
        imax = max(1, min(size(e, 1), i - 1))
        imin = max(1, min(size(e, 1), i - klim))
        j = p2[k+1]
        lUp = min(2 * maxlag + 1, j + 1)
        lDown = max(1, j - 1)
        cUp = sum(e[imin+1:imax, lUp]) + c[imin, lUp]
        cDown = sum(e[imin+1:imax, lDown]) + c[imin, lDown]
        cMiddle = c[imax, j]
        t = findmin([cUp, cDown, cMiddle])
        if t[2] == 1
            p2[imin:k] .= lUp
            p1[imin:k] = imin:k
            k = imin - 1
        elseif t[2] == 2
            p2[imin:k] .= lDown
            p1[imin:k] = imin:k
            k = imin - 1
        else
            p2[k] = j
            p1[k] = k
            k -= 1
        end
    end
    p2 = min.(max.(p2 .+ p1 .- maxlag .- 1, 1), length(p1))
    return (p1, p2)
end # trace

function misfit(p::Setting, m::Vector)
    e = errormap(p["dtw_record"], p["dtw_greenfun"], m, p["dtw_Imaxlag"])
    c = cumulate(e, p["dtw_Iklim"])
    t = c[end, 1]
    for i = 1:size(c, 2)
        if t > c[end, i]
            t = c[end, i]
        end
    end
    return t
end

function detail(p::Setting, m::Vector)
    e = errormap(p["dtw_record"], p["dtw_greenfun"], m, p["dtw_Imaxlag"])
    c = cumulate(e, p["dtw_Iklim"])
    path = trace(e, c, p["dtw_Iklim"])
    return path
end
end # module DTW
