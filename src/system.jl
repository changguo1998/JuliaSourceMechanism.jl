function tomlsortby(e)
    return String(e)
end

"""
parsesac

read sac head variables and return them in an array of Dict
"""
function parsesac(dirname::AbstractString)
    @assert isdir(dirname) "$dirname not exist"
    fs = readdir(dirname)
    metas = Vector{Dict{String,Any}}(undef, length(fs))
    for i in eachindex(fs)
        # h = Seis.readsachead(normpath(dirname, fs[i]))
        h = open(SeisTools.SAC.readhead, normpath(dirname, fs[i]))
        metas[i] = Dict{String,Any}()
        metas[i]["network"] = h["knetwk"]
        metas[i]["station"] = h["kstnm"]
        metas[i]["component"] = String([uppercase(last(h["kcmpnm"]))])
        metas[i]["meta_lon"] = h["stlo"]
        metas[i]["meta_lat"] = h["stla"]
        metas[i]["meta_el"] = h["stel"]
        metas[i]["meta_dt"] = h["delta"]
        metas[i]["meta_btime"] = DateTime(h["nzyear"], 1, 1, h["nzhour"], h["nzmin"], h["nzsec"], h["nzmsec"]) +
                                 Day(h["nzjday"] - 1)
        metas[i]["meta_file"] = fs[i]
    end
    return metas
end

"""
mergemeta

merge related dict in x and y
x and y should be same leng and are Array{Dict, 1}
f is function to pick tag. For example, f(x) = x["network"]*x["station"]
"""
function mergemeta(f::Function, x...)
    if length(x) == 1
        return x[1]
    end
    @assert length(unique(length.(x))) == 1 "meta vectors must be equal"
    tagarray = map(t -> map(f, t), x)
    N = length(tagarray[1])
    T = length(tagarray)
    sortrange = zeros(Int, N, T)
    ref = tagarray[1]
    for i = 1:N, t = 1:T
        sortrange[i, t] = findfirst(tagarray[t] .== ref[i])
    end
    merged = Vector{Dict{String,Any}}(undef, N)
    for n = 1:N
        merged[n] = Dict{String,Any}()
        for t = 2:T
            merge!(merged[n], x[T][sortrange[n, t]])
        end
    end
    return merged
end

"""
mreplace(s::AbstractString, pat...) -> String

replace contents of String `s` with multipul patterns
"""
function mreplace(s::AbstractString, pat...)
    if isempty(pat)
        return s
    else
        t = s
        for p in pat
            t = replace(t, p)
        end
        return t
    end
end

function checkkeyexist!(io::IO, buf::IOBuffer, kk::Tuple, ks, indent::Int, msg::String)
    for k in kk
        if k in ks
            printstyled(io, "    "^indent, k, " OK\n"; color = :green)
        else
            printstyled(io, "    "^indent, k, " ERROR\n"; color = :red)
            println(buf, "$k " * msg)
        end
    end
    return nothing
end

function checkconfiguration!(io::IO, buf::IOBuffer, env::Setting, modules::Vector{Module})
    indent = 1
    if "algorithm" in keys(env)
        printstyled(io, "algorithm OK\n"; color = :green)
        checkkeyexist!(io, buf, ("misfit", "weight", "searchdepth"), keys(env["algorithm"]), indent,
                       "in algorithm is missed.")
    else
        printstyled(io, "algorithm ERROR\n"; color = :red)
        println(buf, "algorithm is missed")
    end
    if "event" in keys(env)
        printstyled(io, "event OK\n"; color = :green)
        checkkeyexist!(io, buf, ("longitude", "latitude", "depth", "magnitude", "origintime"), keys(env["event"]),
                       indent, "in event is missed.")
    else
        printstyled(io, "event ERROR\n"; color = :red)
        println(buf, "event is missed")
    end
    if "stations" in keys(env)
        indent += 1
        if isempty(env["stations"])
            printstyled(io, "stations ERROR\n"; color = :red)
            println(buf, "stations is empty")
        else
            printstyled(io, "stations OK\n"; color = :green)
            KEYS = [["meta_btime", "meta_dt", "meta_el", "meta_file", "meta_lat", "meta_lon", "base_azimuth",
                     "base_distance", "base_trim", "phases"]
                    Green.properties]
            KEYS = Tuple(KEYS)
            MKEYS = ["type", "at", "tt"]
            for m in modules
                append!(MKEYS, m.properties)
            end
            MKEYS = Tuple(MKEYS)
            tbuf2 = IOBuffer()
            for s in env["stations"]
                stag = "    "^(indent - 1) * s["network"] * "." * s["station"] * "." * s["component"]
                checkkeyexist!(devnull, tbuf2, KEYS, keys(s), indent,
                               "in $(s["network"]*"."*s["station"]*"."*s["component"]) is missed")
                if "phases" in keys(s)
                    for p in s["phases"]
                        checkkeyexist!(devnull, tbuf2, MKEYS, keys(p), indent,
                                       "in $(s["network"]*"."*s["station"]*"."*s["component"]).$(p["type"]) is missed")
                    end
                end
                tstr = take!(tbuf2)
                if isempty(tstr)
                    printstyled(io, stag, " OK\n"; color = :green)
                else
                    printstyled(io, stag, " ERROR\n"; color = :red)
                end
                write(buf, tstr)
            end
        end
    else
        printstyled(io, "stations ERROR"; color = :red)
        println(buf, "stations is missed")
    end
    return nothing
end

function checkconfiguration(io::IO, env::Setting, modules::Vector{Module})
    buf = IOBuffer()
    checkconfiguration!(io, buf, env, modules)
    return String(take!(buf))
end

function checkconfiguration(env::Setting, modules::Vector{Module})
    buf = IOBuffer()
    checkconfiguration!(stdout, buf, env, modules)
    return String(take!(buf))
end

function buildstationconfiguration(root::AbstractString, event::Setting)
    sacmeta = parsesac(normpath(root, "sac"))
    for s in sacmeta
        (dist, az, _) = SeisTools.Geodesy.distance(event["latitude"], event["longitude"], s["meta_lat"], s["meta_lon"])
        s["base_distance"] = dist
        s["base_azimuth"] = az
        s["base_trim"] = [s["meta_btime"], s["meta_btime"] + Minute(10)]
        s["phases"] = Dict{String,Any}[]
    end
    return sacmeta
end

function init_event(dir::AbstractString)
    if !isdir(dir)
        return nothing
    end
    for f in filter(v -> endswith(v, ".jl") || endswith(v, ".m") || endswith(v, ".fig"),
                    readdir(joinpath(@__DIR__, "..", "example")))
        cp(joinpath(@__DIR__, "..", "example", f), joinpath(dir, f); force = true)
    end
    return nothing
end
