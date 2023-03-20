module RelShift
using Dates, DSP, Statistics, LinearAlgebra
import JuliaSourceMechanism: Setting

tags = ("relshift", "rshift", "RelShift", "Rshift", "rsft", "Rsft")
properties = ["relshift_dt", "relshift_order", "relshift_band", "relshift_maxlag", "relshift_trim"]

function weight(p::Setting, s::Setting, e::Setting)
    if "relshift_weight" in keys(p)
        return p["relshift_weight"]
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
    if station["component"] != "Z"
        return nothing
    end
    cmpsloc = findall(x -> x["network"] == station["network"] && x["station"] == station["station"], env["stations"])
    cmplist = String[]
    for sta in env["stations"][cmpsloc]
        c = uppercase(String(sta["component"]))
        p_idx = findfirst(p -> p["type"] == phase["type"], sta["phases"])
        p = sta["phases"][p_idx]
        w = deepcopy(sta["base_record"])
        w_resample = resample(w, sta["meta_dt"] / phase["relshift_dt"])
        fltr = digitalfilter(Bandpass(phase["relshift_band"][1], phase["relshift_band"][2];
                                      fs = 1 / phase["relshift_dt"]),
                             Butterworth(phase["relshift_order"]))
        w_filt = filtfilt(fltr, w_resample)
        (_, w_trim, _) = cut(w_filt, sta["base_begintime"],
                             phase["at"] + Millisecond(round(Int, phase["relshift_trim"][1] * 1e3)),
                             phase["at"] + Millisecond(round(Int, phase["relshift_trim"][2] * 1e3)),
                             Millisecond(round(Int, 1e3 * phase["relshift_dt"])))
        nm = norm(w_trim)
        w_trim ./= nm
        g = deepcopy(sta["green_fun"])
        g_resample = resample(g, sta["green_dt"] / phase["relshift_dt"]; dims = 1)
        g_filt = filtfilt(fltr, g_resample)
        (_, g_trim, _) = cut(g_filt, station["base_begintime"],
                             station["base_begintime"] +
                             Millisecond(round(Int, (phase["tt"] + phase["relshift_trim"][1]) * 1e3)),
                             station["base_begintime"] +
                             Millisecond(round(Int, (phase["tt"] + phase["relshift_trim"][2]) * 1e3)),
                             Millisecond(round(Int, 1e3 * phase["relshift_dt"])); fillval = 0.0)
        txcorr = permutedims(xcorr(w_trim, g_trim;
                                   maxlag = round(Int, phase["relshift_maxlag"] / phase["relshift_dt"])))
        phase["relshift_relation_"*c] = txcorr
        amp = zeros(6, 6)
        for i = 1:6, j = 1:6, k in axes(g_trim, 1)
            amp[i, j] += g_trim[k, i] * g_trim[k, j]
        end
        phase["relshift_record_"*c] = w_trim
        phase["relshift_greenfun_"*c] = g_trim
        phase["relshift_shiftbase_"*c] = size(g_trim, 1)
        push!(cmplist, c)
    end
    phase["relshift_componentlist"] = Tuple(cmplist)
    return nothing
end

function misfit(p::Setting, m::Vector)
    maxv = 0.0
    lcmp = length(p["relshift_componentlist"])
    idx = zeros(Int, lcmp)
    for c = 1:lcmp
        tr = p["relshift_relation_"*p["relshift_componentlist"][c]]
        idx[c] = 0
        for j in axes(tr, 2)
            local tv = 0.0
            for i = 1:6
                tv += tr[i, j] * m[i]
            end
            if tv > maxv
                maxv = tv
                idx[c] = j
            end
        end
    end
    var = 0.0
    for c = 1:lcmp-1
        var += abs2(idx[c] - idx[c+1])
    end
    var += abs2(idx[end] - idx[1])
    var /= lcmp
    return sqrt(var) * p["relshift_dt"] * 2.0 / sqrt(3.0) / p["relshift_maxlag"]
end
end
