using JuliaSourceMechanism, JLD2, LinearAlgebra, Printf

@info "Load buffer"
(env, status) = let
    t = load(joinpath(@__DIR__, "..", "data.jld2"))
    (t["env"], t["status"])
end

@info "Copy data"
nenv = Setting()
nenv["dataroot"] = abspath("..")
nenv["algorithm"] = env["algorithm"]
nenv["event"] = env["event"]
nenv["stations"] = Setting[]
for s in env["stations"]
    if s["network"] * "." * s["station"] ∈ status["stationlist_selected"]
        t = deepcopy(s)
        push!(nenv["stations"], t)
    end
end

# @info "Change parameter"
# nenv["algorithm"]["misfit"] = ["xcorr", "pol"]
# nenv["algorithm"]["weight"] = [1.0, 0.4]
# for s in nenv["stations"]
#     s["green_tsource"] = 0.04
#     s["green_t0"] = -2.0 * s["green_tsource"]
#     for p in s["phases"]
#         dt = s["meta_dt"]
#         tl = p["xcorr_trim"][2] - p["xcorr_trim"][1]
#         while (dt + 1e-3) * 200 < tl
#             dt += 1e-3
#         end
#         p["xcorr_dt"] = dt
#     end
# end

misfits = Module[]
for m in nenv["algorithm"]["misfit"], f in [XCorr, Polarity, PSR, DTW, AbsShift, RelShift]
    if m in f.tags
        push!(misfits, f)
    end
end

@info "Run"
preprocess!(nenv, misfits; warn = false)
(sdr, phaselist, misfit, misfitdetail) = inverse!(nenv, misfits, Grid)
weight = normalize(map(x -> x[1].weight(x[2], env, env), phaselist), 1)
totalmisfit = replace(misfitdetail, NaN => 0.0) * weight
(minval, minidx) = findmin(totalmisfit)
mech = sdr[minidx]

result = Setting()
for s in nenv["stations"]
    sta = s["network"] * "." * s["station"]
    if !(sta ∈ keys(result))
        result[sta] = Setting()
    end
    c = s["component"]
    result[sta]["dist"] = s["base_distance"]
    result[sta][c] = Setting()
    for p in s["phases"]
        t = Setting()
        t["xcorr_rec"] = normalize(p["xcorr_record"], Inf)
        t["xcorr_syn"] = normalize(p["xcorr_greenfun"] * dc2ts(mech), Inf)
        t["xcorr_shift"] = XCorr.detail(p, dc2ts(mech))
        t["xcorr_dt"] = p["xcorr_dt"]
        t["polarity_rec"] = p["polarity_obs"]
        if p["type"] == "P"
            t["polarity_syn"] = sign(sum(p["polarity_syn"] .* dc2ts(mech)))
        end
        result[sta][c][p["type"]] = t
    end
end

@info "Save"
jldsave(joinpath(@__DIR__, "result.jld2"); mech=mech, result=result)

@info "Plot"
include("plot.jl")
