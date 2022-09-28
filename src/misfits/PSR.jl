module PSR

using Dates, Statistics
import JuliaSourceMechanism: Setting, trim

tags = ("psr", "PSR")
properties = ["psr_trimp", "psr_trims"]

function weight(p::Setting, s::Setting, e::Setting)
    if "psr_weight" ∈ keys(p)
        return p["psr_weight"]
    else
        idx = findfirst(x -> x ∈ tags, e["algorithm"]["misfit"])
        if isnothing(idx)
            return 0.0
        else
            return e["algorithm"]["weight"][idx]
        end
    end
end

function preprocess!(phase::Setting, station::Setting, env::Setting)
    @debug "prepare $(tags[1]) data for station: $(station["network"]).$(station["station"]).$(station["component"]) phase: $(phase["type"])"
    if phase["type"] == "S"
        phase["psr_obs"] = NaN
        return nothing
    else
        typelist = map(x -> x["type"], station["phases"])
        if length(typelist) < 2
            phase["psr_obs"] = NaN
            return nothing
        end
        idx = findfirst(typelist .== "S")
        w = deepcopy(station["base_record"])
        wp = trim(w, station["base_begintime"], phase["at"] + Millisecond(round(Int, phase["psr_trimp"][1])),
                  phase["at"] + Millisecond(round(Int, phase["psr_trimp"][2])), station["meta_dt"])
        ws = trim(w, station["base_begintime"],
                  station["phases"][idx]["at"] + Millisecond(round(Int, phase["psr_trims"][1])),
                  station["phases"][idx]["at"] + Millisecond(round(Int, phase["psr_trims"][2])), station["meta_dt"])
        phase["psr_obs"] = 10.0 * log10(mean(abs2, ws) / mean(abs2, wp))
        g = deepcopy(station["green_fun"])
        o = env["event"]["origintime"]
        g_trimp = trim(g, o, o + Millisecond(round(Int, (phase["tt"] + phase["psr_trimp"][1]) * 1e3)),
                       o + Millisecond(round(Int, (phase["tt"] + phase["psr_trimp"][2]) * 1e3)), station["green_dt"])
        g_trims = trim(g, o, o + Millisecond(round(Int, (phase["tt"] + phase["psr_trims"][1]) * 1e3)),
                       o + Millisecond(round(Int, (phase["tt"] + phase["psr_trims"][2]) * 1e3)), station["green_dt"])
        amp_p = zeros(6, 6)
        amp_s = zeros(6, 6)
        for i = 1:6, j = 1:6
            for k = 1:size(g_trimp, 1)
                amp_p[i, j] += g_trimp[k, i] * g_trimp[k, j]
            end
            for k = 1:size(g_trims, 1)
                amp_s[i, j] += g_trims[k, i] * g_trims[k, j]
            end
        end
        amp_s ./= size(g_trimp, 1) / size(g_trims, 1)
        phase["psr_ampP"] = amp_p
        phase["psr_ampS"] = amp_s
        return nothing
    end
end

function misfit(phase::Setting, m::Vector)
    if isnan(phase["psr_obs"])
        return NaN
    else
        ap = 0.0
        for i = 1:6, j = 1:6
            ap += phase["psr_ampP"][i, j] * m[i] * m[j]
        end
        as = 0.0
        for i = 1:6, j = 1:6
            as += phase["psr_ampS"][i, j] * m[i] * m[j]
        end
        return abs(10 * log10(as / ap) - phase["psr_obs"])
    end
end
end
