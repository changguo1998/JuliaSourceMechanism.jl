module AbsShift
using Dates, DSP, Statistics, LinearAlgebra, SeisTools.DataProcess
import JuliaSourceMechanism: Setting

tags = ("absshift", "ashift", "AbsShift", "Ashift", "asft", "Asft")
properties = ["absshift_dt", "absshift_order", "absshift_band", "absshift_maxlag", "absshift_trim"]

function weight(p::Setting, s::Setting, e::Setting)
    if "absshift_weight" in keys(p)
        return p["absshift_weight"]
    else
        idx = findfirst(t -> t in tags, e["algorithm"]["misfit"])
        return e["algorithm"]["weight"][idx]
    end
end

"""
"""
function xcorr(u::VecOrMat, v::VecOrMat; maxlag::Union{Int,Nothing} = nothing)
    Lu = size(u, 1)
    Lv = size(v, 1)
    Wu = size(u, 2)
    Wv = size(v, 2)
    if isnothing(maxlag)
        maxlag = min(Lu, Lv) - 1
    end
    r = zeros(2 * maxlag + 1, Wu * Wv)
    for i in axes(r, 1), j in axes(r, 2)
        s = i - maxlag - 1
        (ju, jv) = divrem(j - 1, Wv)
        ju += 1
        jv += 1
        minv = max(1, 1 - s)
        maxv = min(Lv, Lu - s)
        minu = max(1, minv + s)
        r[i, j] = 0.0
        for l = 0:(maxv-minv)
            r[i, j] += u[minu+l, ju] * v[minv+l, jv]
        end
    end
    return r
end

function preprocess!(phase::Setting, station::Setting, env::Setting)
    w = deepcopy(station["base_record"])
    w_resample = resample(w, station["meta_dt"] / phase["absshift_dt"])
    fltr = digitalfilter(Bandpass(phase["absshift_band"][1], phase["absshift_band"][2]; fs = 1 / phase["absshift_dt"]),
                         Butterworth(phase["absshift_order"]))
    w_filt = filtfilt(fltr, w_resample)
    # w_trim = trim(w_filt, station["base_begintime"],
    #               phase["at"] + Millisecond(round(Int, phase["absshift_trim"][1] * 1e3)),
    #               phase["at"] + Millisecond(round(Int, phase["absshift_trim"][2] * 1e3)), phase["absshift_dt"])
    (_, w_trim, _) = cut(w_filt, station["base_begintime"],
                         phase["at"] + Millisecond(round(Int, phase["absshift_trim"][1] * 1e3)),
                         phase["at"] + Millisecond(round(Int, phase["absshift_trim"][2] * 1e3)),
                         Millisecond(round(Int, phase["absshift_dt"] * 1e3)))
    nm = norm(w_trim)
    w_trim ./= nm
    g = deepcopy(station["green_fun"])
    g_resample = resample(g, station["green_dt"] / phase["absshift_dt"]; dims = 1)
    g_filt = filtfilt(fltr, g_resample)
    (_, g_trim, _) = cut(g_filt, station["base_begintime"],
                         station["base_begintime"] +
                         Millisecond(round(Int, (phase["tt"] + phase["absshift_trim"][1]) * 1e3)),
                         station["base_begintime"] +
                         Millisecond(round(Int, (phase["tt"] + phase["absshift_trim"][2]) * 1e3)),
                         Millisecond(round(Int, phase["absshift_dt"] * 1e3)); fillval = 0.0)
    txcorr = permutedims(xcorr(w_trim, g_trim; maxlag = round(Int, phase["absshift_maxlag"] / phase["absshift_dt"])))
    phase["absshift_relation"] = txcorr
    amp = zeros(6, 6)
    for i = 1:6, j = 1:6, k in axes(g_trim, 1)
        amp[i, j] += g_trim[k, i] * g_trim[k, j]
    end
    phase["absshift_record"] = w_trim
    phase["absshift_greenfun"] = g_trim
    phase["absshift_shiftbase"] = size(g_trim, 1)
    return nothing
end

function misfit(p::Setting, m::Vector)
    maxv = 0.0
    tr = deepcopy(p["absshift_relation"])
    idx = 0
    for j in axes(tr, 2)
        local tv = 0.0
        for i = 1:6
            tv += tr[i, j] * m[i]
        end
        if tv > maxv
            maxv = tv
            idx = j
        end
    end
    return abs2((idx - p["absshift_shiftbase"]) * p["absshift_dt"]) / ["absshift_maxlag"]
end
end
