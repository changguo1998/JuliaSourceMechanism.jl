module XCorr
using Dates, DSP, Statistics, LinearAlgebra, SeisTools.DataProcess
import JuliaSourceMechanism: Setting

tags = ("XCorr", "xcorr")
properties = ["xcorr_dt", "xcorr_order", "xcorr_band", "xcorr_maxlag", "xcorr_trim"]

function _Second(t::Real, prec::Type=Millisecond)
    return prec(round(Int, t*(prec(Second(1))/prec(1))))
end

function weight(p::Setting, s::Setting, e::Setting)
    if "xcorr_weight" in keys(p)
        return p["xcorr_weight"]
    else
        idx = findfirst(t -> t in tags, e["algorithm"]["misfit"])
        return e["algorithm"]["weight"][idx]
    end
end

"""
"""
function _xcorr(u::VecOrMat, v::VecOrMat; maxlag::Union{Int,Nothing} = nothing)
    Lu = size(u, 1)
    Lv = size(v, 1)
    Wu = size(u, 2)
    Wv = size(v, 2)
    if isnothing(maxlag)
        maxlag = min(Lu, Lv) - 1
    end
    if maxlag > Lu + Lv - 1
        @warn "maxlag $(maxlag) is too large, set to $(Lu + Lv - 1)"
        maxlag = Lu + Lv - 1
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
        for l = 0:(maxv-minv)
            r[i, j] += u[minu+l, ju] * v[minv+l, jv]
        end
    end
    return r
end

function preprocess!(phase::Setting, station::Setting, env::Setting)
    @debug "prepare $(tags[1]) data for station: $(station["network"]).$(station["station"]).$(station["component"]) \
    phase: $(phase["type"])"
    w = deepcopy(station["base_record"])
    w_resample = DataProcess.resample(w, station["meta_dt"] / phase["xcorr_dt"])
    fltr = digitalfilter(Bandpass(phase["xcorr_band"][1], phase["xcorr_band"][2]),
                         Butterworth(phase["xcorr_order"]); fs = 1 / phase["xcorr_dt"])
    w_filt = filtfilt(fltr, w_resample)
    (_, w_trim, _) = cut(w_filt, station["base_begintime"],
                         phase["at"] + _Second(phase["xcorr_trim"][1]),
                         phase["at"] + _Second(phase["xcorr_trim"][2]),
                         _Second(phase["xcorr_dt"]))
    nm = norm(w_trim)
    w_trim ./= nm
    g = deepcopy(station["green_fun"])
    g_resample = DataProcess.resample(g, station["green_dt"] / phase["xcorr_dt"])
    g_filt = filtfilt(fltr, g_resample)
    g_trim = zeros(length(w_trim), 6)
    _gshift = _Second(phase["tt"]+phase["xcorr_trim"][1]) / _Second(phase["xcorr_dt"]) |> _t->round(Int, _t)
    cut!(g_trim, g_filt, _gshift, 1, length(w_trim))
    # (_, g_trim, _) = cut(g_filt, station["base_begintime"],
    #                      station["base_begintime"] +
    #                      Millisecond(round(Int, (phase["tt"] + phase["xcorr_trim"][1]) * 1e3)),
    #                      station["base_begintime"] +
    #                      Millisecond(round(Int, (phase["tt"] + phase["xcorr_trim"][2]) * 1e3)),
    #                      Millisecond(round(Int, 1e3 * phase["xcorr_dt"])); fillval = 0.0)
    txcorr = permutedims(_xcorr(w_trim, g_trim; maxlag = round(Int, phase["xcorr_maxlag"] / phase["xcorr_dt"])))
    phase["xcorr_relation"] = txcorr
    amp = zeros(6, 6)
    for i = 1:6, j = 1:6, k in axes(g_trim, 1)
        amp[i, j] += g_trim[k, i] * g_trim[k, j]
    end
    phase["xcorr_synamp"] = amp
    phase["xcorr_record"] = w_trim
    phase["xcorr_greenfun"] = g_trim
    return nothing
end

function misfit(p::Setting, m::Vector)
    amp = 0.0
    for j = 1:6, i = 1:6
        amp += p["xcorr_synamp"][i, j] * m[i] * m[j]
    end
    amp = sqrt(amp)
    maxv = 0.0
    tv = 0.0
    tr = deepcopy(p["xcorr_relation"])
    for j in axes(tr, 2)
        tv = 0.0
        for i = 1:6
            tv += tr[i, j] * m[i]
        end
        if tv > maxv
            maxv = tv
        end
    end
    return sqrt(1.0 - maxv / amp)
end

function detail(p::Setting, m::Vector)
    amp = 0.0
    for j = 1:6, i = 1:6
        amp += p["xcorr_synamp"][i, j] * m[i] * m[j]
    end
    amp = sqrt(amp)
    maxv = 0.0
    tv = 0.0
    iv = 0
    tr = deepcopy(p["xcorr_relation"])
    for j in axes(tr, 2)
        tv = 0.0
        for i = 1:6
            tv += tr[i, j] * m[i]
        end
        if tv > maxv
            maxv = tv
            iv = j
        end
    end
    n = size(p["xcorr_relation"], 2)
    n = round(Int, (n + 1) / 2)
    return (iv - n) * p["xcorr_dt"]
end

end # module JINV_xcorr
