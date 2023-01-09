module CAP
using Dates, Statistics, LinearAlgebra, SeisTools.DataProcess
import JuliaSourceMechanism: Setting

tags = ("cap", "CAP")
properties = ["cap_dt", "cap_order", "cap_band_pnl", "cap_band_all", "cap_r0",
              "cap_p_pnl", "cap_p_rayleigh", "cap_p_love",
              "cap_pnl_maxlag", "cap_all_maxlag", "cap_pnl_trim", "cap_all_trim"]

function weight(p::Setting, s::Setting, e::Setting)
    if "cap_weight" ∈ keys(p)
        return p["cap_weight"]
    else
        idx = findfirst(x -> x ∈ tags, e["algorithm"]["misfit"])
        if isnothing(idx)
            return 0.0
        else
            return e["algorithm"]["weight"][idx]
        end
    end
end

function skip(p::Setting)
    return p["cap_skip"]
end

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
    if (station["component"] != "Z") || (phase["type"] != "P")
        phase["cap_skip"] = true
        return nothing
    end
    phase["cap_skip"] = false
    cmp_e = findfirst(x -> x["network"] == station["network"] &&
                             x["station"] == station["station"] &&
                             x["component"] == "E", env["stations"])
    cmp_n = findfirst(x -> x["network"] == station["network"] &&
                             x["station"] == station["station"] &&
                             x["component"] == "N", env["stations"])
    we = deepcopy(env["stations"][cmp_e]["base_record"])
    wn = deepcopy(env["stations"][cmp_n]["base_record"])
    wz = deepcopy(station["base_record"])
    T = [sind(station["base_azimuth"]) -cosd(station["base_azimuth"]);
         cosd(station["base_azimuth"]) sind(station["base_azimuth"])]
    w = [we wn] * T
    wr = w[:, 1]
    wt = w[:, 2]

    wr_resample = resample(wr, station["meta_dt"] / phase["cap_dt"])
    wt_resample = resample(wt, station["meta_dt"] / phase["cap_dt"])
    wz_resample = resample(wz, station["meta_dt"] / phase["cap_dt"])
    # fltr = digitalfilter(Bandpass(phase["xcorr_band"][1], phase["xcorr_band"][2]; fs = 1 / phase["xcorr_dt"]),
    #                      Butterworth(phase["xcorr_order"]))
    # we_filt = filtfilt(fltr, we_resample)
    # wn_filt = filtfilt(fltr, wn_resample)
    # wz_filt = filtfilt(fltr, wz_resample)
    wr_filt_pnl = bandpass(wr_resample, phase["cap_band_pnl"][1], phase["cap_band_pnl"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    wz_filt_pnl = bandpass(wz_resample, phase["cap_band_pnl"][1], phase["cap_band_pnl"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    wr_filt_all = bandpass(wr_resample, phase["cap_band_all"][1], phase["cap_band_all"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    wt_filt_all = bandpass(wt_resample, phase["cap_band_all"][1], phase["cap_band_all"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    wz_filt_all = bandpass(wz_resample, phase["cap_band_all"][1], phase["cap_band_all"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    (_, wr_pnl_trim, _) = cut(wr_filt_pnl, station["base_begintime"],
                              phase["at"] + Millisecond(round(Int, phase["cap_pnl_trim"][1] * 1e3)),
                              phase["at"] + Millisecond(round(Int, phase["cap_pnl_trim"][2] * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))
    (_, wz_pnl_trim, _) = cut(wz_filt_pnl, station["base_begintime"],
                              phase["at"] + Millisecond(round(Int, phase["cap_pnl_trim"][1] * 1e3)),
                              phase["at"] + Millisecond(round(Int, phase["cap_pnl_trim"][2] * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))

    (_, wr_all_trim, _) = cut(wr_filt_all, station["base_begintime"],
                              phase["at"] + Millisecond(round(Int, phase["cap_all_trim"][1] * 1e3)),
                              phase["at"] + Millisecond(round(Int, phase["cap_all_trim"][2] * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))
    (_, wt_all_trim, _) = cut(wt_filt_all, station["base_begintime"],
                              phase["at"] + Millisecond(round(Int, phase["cap_all_trim"][1] * 1e3)),
                              phase["at"] + Millisecond(round(Int, phase["cap_all_trim"][2] * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))
    (_, wz_all_trim, _) = cut(wz_filt_all, station["base_begintime"],
                              phase["at"] + Millisecond(round(Int, phase["cap_all_trim"][1] * 1e3)),
                              phase["at"] + Millisecond(round(Int, phase["cap_all_trim"][2] * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))
    # nm = norm(w_trim)
    # w_trim ./= nm
    ge = deepcopy(env["stations"][cmp_e]["green_fun"])
    gn = deepcopy(env["stations"][cmp_n]["green_fun"])
    gz = deepcopy(station["green_fun"])
    gr = ge .* sind(station["base_azimuth"]) .+ gn * cosd(station["base_azimuth"])
    gt = ge .* -cosd(station["base_azimuth"]) .+ gn * sind(station["base_azimuth"])
    gr_resample = resample(gr, station["green_dt"] / phase["cap_dt"])
    gt_resample = resample(gt, station["green_dt"] / phase["cap_dt"])
    gz_resample = resample(gz, station["green_dt"] / phase["cap_dt"])
    # g_filt = filtfilt(fltr, g_resample)
    gr_filt_pnl = bandpass(gr_resample, phase["cap_band_pnl"][1], phase["cap_band_pnl"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    gz_filt_pnl = bandpass(gz_resample, phase["cap_band_pnl"][1], phase["cap_band_pnl"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    gr_filt_all = bandpass(gr_resample, phase["cap_band_all"][1], phase["cap_band_all"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    gt_filt_all = bandpass(gt_resample, phase["cap_band_all"][1], phase["cap_band_all"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    gz_filt_all = bandpass(gz_resample, phase["cap_band_all"][1], phase["cap_band_all"][2],
                           1.0 / phase["cap_dt"]; n = phase["cap_order"])
    (_, gr_pnl_trim, _) = cut(gr_filt_pnl, station["base_begintime"],
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_pnl_trim"][1]) * 1e3)),
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_pnl_trim"][2]) * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))
    (_, gz_pnl_trim, _) = cut(gz_filt_pnl, station["base_begintime"],
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_pnl_trim"][1]) * 1e3)),
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_pnl_trim"][2]) * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))

    (_, gr_all_trim, _) = cut(gr_filt_all, station["base_begintime"],
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_all_trim"][1]) * 1e3)),
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_all_trim"][2]) * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))
    (_, gt_all_trim, _) = cut(gt_filt_all, station["base_begintime"],
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_all_trim"][1]) * 1e3)),
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_all_trim"][2]) * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))
    (_, gz_all_trim, _) = cut(gz_filt_all, station["base_begintime"],
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_all_trim"][1]) * 1e3)),
                              station["base_begintime"] +
                              Millisecond(round(Int, (phase["tt"] + phase["cap_all_trim"][2]) * 1e3)),
                              Millisecond(round(Int, 1e3 * phase["cap_dt"])))
    # txcorr = permutedims(_xcorr(w_trim, g_trim; maxlag = round(Int, phase["xcorr_maxlag"] / phase["xcorr_dt"])))
    txc_r_pnl = _xcorr(wr_pnl_trim, gr_pnl_trim; maxlag = round(Int, phase["cap_pnl_maxlag"] / phase["cap_dt"])) |>
                permutedims
    txc_z_pnl = _xcorr(wz_pnl_trim, gz_pnl_trim; maxlag = round(Int, phase["cap_pnl_maxlag"] / phase["cap_dt"])) |>
                permutedims
    txc_r_all = _xcorr(wr_all_trim, gr_all_trim; maxlag = round(Int, phase["cap_all_maxlag"] / phase["cap_dt"])) |>
                permutedims
    txc_t_all = _xcorr(wt_all_trim, gt_all_trim; maxlag = round(Int, phase["cap_all_maxlag"] / phase["cap_dt"])) |>
                permutedims
    txc_z_all = _xcorr(wz_all_trim, gz_all_trim; maxlag = round(Int, phase["cap_all_maxlag"] / phase["cap_dt"])) |>
                permutedims
    # phase["xcorr_relation"] = txcorr
    phase["cap_relation_r_pnl"] = txc_r_pnl
    phase["cap_relation_z_pnl"] = txc_z_pnl
    phase["cap_relation_r_all"] = txc_r_all
    phase["cap_relation_t_all"] = txc_t_all
    phase["cap_relation_z_all"] = txc_z_all
    amp_r_pnl = zeros(6, 6)
    amp_z_pnl = zeros(6, 6)
    amp_r_all = zeros(6, 6)
    amp_t_all = zeros(6, 6)
    amp_z_all = zeros(6, 6)
    for i = 1:6, j = 1:6
        for k in axes(gr_pnl_trim, 1)
            amp_r_pnl[i, j] += gr_pnl_trim[k, i] * gr_pnl_trim[k, j]
            amp_z_pnl[i, j] += gz_pnl_trim[k, i] * gz_pnl_trim[k, j]
        end
        for k in axes(gr_all_trim, 1)
            amp_r_all[i, j] += gr_all_trim[k, i] * gr_all_trim[k, j]
            amp_t_all[i, j] += gt_all_trim[k, i] * gt_all_trim[k, j]
            amp_z_all[i, j] += gz_all_trim[k, i] * gz_all_trim[k, j]
        end
    end
    # phase["xcorr_synamp"] = amp
    phase["cap_synamp_r_pnl"] = amp_r_pnl
    phase["cap_synamp_z_pnl"] = amp_z_pnl
    phase["cap_synamp_r_all"] = amp_r_all
    phase["cap_synamp_t_all"] = amp_t_all
    phase["cap_synamp_z_all"] = amp_z_all
    phase["cap_record_r_pnl"] = wr_pnl_trim
    phase["cap_record_z_pnl"] = wz_pnl_trim
    phase["cap_record_r_all"] = wr_all_trim
    phase["cap_record_t_all"] = wt_all_trim
    phase["cap_record_z_all"] = wz_all_trim
    phase["cap_f2_r_pnl"] = sum(abs2, wr_pnl_trim)
    phase["cap_f2_z_pnl"] = sum(abs2, wz_pnl_trim)
    phase["cap_f2_r_all"] = sum(abs2, wr_all_trim)
    phase["cap_f2_t_all"] = sum(abs2, wt_all_trim)
    phase["cap_f2_z_all"] = sum(abs2, wz_all_trim)
    # phase["xcorr_greenfun"] = g_trim
    phase["cap_greenfun_r_pnl"] = gr_pnl_trim
    phase["cap_greenfun_z_pnl"] = gz_pnl_trim
    phase["cap_greenfun_r_all"] = gr_all_trim
    phase["cap_greenfun_t_all"] = gt_all_trim
    phase["cap_greenfun_z_all"] = gz_all_trim
    phase["cap_distance"] = station["base_distance"]
    return nothing
end

function rec2(phase::Setting)
    return phase["cap_f2_r_pnl"] +
           phase["cap_f2_z_pnl"] +
           phase["cap_f2_r_all"] +
           phase["cap_f2_t_all"] +
           phase["cap_f2_z_all"]
end

function _amp2(mat::Matrix, m::Vector)
    v = 0.0
    for i = 1:6, j = 1:6
        v += mat[i, j] * m[i] * m[j]
    end
    return v
end

function syn2(phase::Setting, m::Vector)
    map = phase["cap_synamp_r_pnl"] +
          phase["cap_synamp_z_pnl"] +
          phase["cap_synamp_r_all"] +
          phase["cap_synamp_t_all"] +
          phase["cap_synamp_z_all"]
    return _amp2(map, m)
end

function _maxcorr(xcor::Matrix, m::Vector)
    maxv = -Inf
    maxi = 0
    for j in axes(xcor, 2)
        v = 0.0
        for i = 1:6
            v += xcor[i, j] * m[i]
        end
        if v > maxv
            maxv = v
            maxi = j
        end
    end
    return (maxv, maxi)
end

function misfit(p::Setting, m::Vector)
    if p["cap_skip"]
        return 0.0
    end
    (v_r_pnl, _) = _maxcorr(p["cap_relation_r_pnl"], m)
    (v_z_pnl, _) = _maxcorr(p["cap_relation_z_pnl"], m)
    (v_r_all, _) = _maxcorr(p["cap_relation_r_all"], m)
    (v_t_all, _) = _maxcorr(p["cap_relation_t_all"], m)
    (v_z_all, _) = _maxcorr(p["cap_relation_z_all"], m)
    mis_r_pnl = sqrt(p["cap_f2_r_pnl"] + _amp2(p["cap_synamp_r_pnl"], m) - 2 * v_r_pnl)
    mis_z_pnl = sqrt(p["cap_f2_z_pnl"] + _amp2(p["cap_synamp_z_pnl"], m) - 2 * v_z_pnl)
    mis_r_all = sqrt(p["cap_f2_r_all"] + _amp2(p["cap_synamp_r_all"], m) - 2 * v_r_all)
    mis_t_all = sqrt(p["cap_f2_t_all"] + _amp2(p["cap_synamp_t_all"], m) - 2 * v_t_all)
    mis_z_all = sqrt(p["cap_f2_z_all"] + _amp2(p["cap_synamp_z_all"], m) - 2 * v_z_all)
    w_pnl = (p["cap_distance"] / p["cap_r0"])^p["cap_p_pnl"]
    w_psv = (p["cap_distance"] / p["cap_r0"])^p["cap_p_rayleigh"]
    w_sh = (p["cap_distance"] / p["cap_r0"])^p["cap_p_love"]
    return w_pnl * (mis_r_pnl + mis_z_pnl) + w_psv * (mis_r_all + mis_z_all) + w_sh * mis_t_all
end

function detail(p::Setting, m::Vector)
    (_, i_r_pnl) = _maxcorr(p["cap_relation_r_pnl"], m)
    (_, i_z_pnl) = _maxcorr(p["cap_relation_z_pnl"], m)
    (_, i_r_all) = _maxcorr(p["cap_relation_r_all"], m)
    (_, i_t_all) = _maxcorr(p["cap_relation_t_all"], m)
    (_, i_z_all) = _maxcorr(p["cap_relation_z_all"], m)
    ref_pnl = round(Int, (size(p["cap_relation_r_pnl"], 2) + 1) / 2)
    ref_all = round(Int, (size(p["cap_relation_r_all"], 2) + 1) / 2)
    return (i_r_pnl - ref_pnl, i_z_pnl - ref_pnl, i_r_all - ref_all, i_t_all - ref_all, i_z_all - ref_all) .*
           p["cap_dt"]
end

end
