module Polarity

using Dates, SeisTools.DataProcess
import JuliaSourceMechanism: Setting

tags = ("polarity", "pol", "Pol", "Polarity")
properties = ["polarity_trim", "polarity_obs"]

function weight(p::Setting, s::Setting, e::Setting)
    if "polarity_weight" ∈ keys(p)
        return p["polarity_weight"]
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
    if phase["type"] != "P"
        return nothing
    end
    g = deepcopy(station["green_fun"])
    o = station["base_begintime"]
    (_, g_trim, _) = cut(g, o, o + Millisecond(round(Int, (phase["tt"] + phase["polarity_trim"][1]) * 1e3)),
                  o + Millisecond(round(Int, (phase["tt"] + phase["polarity_trim"][2]) * 1e3)), 
                  Millisecond(round(Int, station["green_dt"]*1e3)))
    tvec = zeros(6)
    for j = 1:6
        for i = axes(g_trim, 1)
            tvec[j] += g_trim[i, j]
        end
    end
    phase["polarity_syn"] = Tuple(tvec)
    phase["polarity_greenfun"] = g_trim
    return nothing
end

function msign(x::Float64)
    if x > 0.0
        return 1.0
    elseif x < 0.0
        return -1.0
    else
        return 0.0
    end
end

function misfit(phase::Setting, m::Vector)
    if isnan(phase["polarity_obs"])
        return NaN
    end
    t = 0.0
    for i = 1:6
        t += phase["polarity_syn"][i] * m[i]
    end
    t = msign(t)
    if t == phase["polarity_obs"] || iszero(phase["polarity_obs"])
        return 0.0
    else
        return 1.0
    end
end
end # module JINV_polarity
