module InteractiveTest

using Dates, DelimitedFiles, DSP, FFTW, JLD2, MAT, Printf, Statistics, TOML, SeisTools
import JuliaSourceMechanism: Setting, Green, dc2ts, trim, detrendandtaper!, checkconfiguration
include("ArgParser.jl")

function fresh!(env::Setting, misfitmodules::Vector{Module}; warn::Bool = false)
    t = checkconfiguration(devnull, env, misfitmodules)
    if !isempty(t)
        @error "Something wrong while checking. See $(abspath(normpath(pwd(), "check.log"))) for more detail."
        open("check.log", "w") do io
            print(io, t)
        end
        return nothing
    end
    for s in env["stations"]
        (dist, az, _) = SeisTools.Geodesy.distance(env["event"]["latitude"], env["event"]["longitude"], s["meta_lat"], s["meta_lon"])
        s["base_distance"] = dist
        s["base_azimuth"] = az
    end

    @info "Green function"
    if !("dataroot" ∈ keys(env))
        @warn "dataroot not specified"
        return nothing
    end
    # Threads.@threads
    for s in collect(filter(x -> x["component"] == "Z", env["stations"]))
        Green.calculategreenfun(s, env)
    end

    for s in env["stations"]
        # t = Seis.readsac(normpath(env["dataroot"], "sac", s["meta_file"]))
        t = SeisTools.SAC.read(normpath(env["dataroot"], "sac", s["meta_file"]))
        # if s["base_trim"][1] < s["meta_btime"]
        #     trim_bt = s["meta_btime"]
        #     if warn
        #         printstyled("Begin time: ", s["base_trim"][1],
        #                     " of $(s["network"]).$(s["station"]).$(s["component"]) is too early, set to: ", trim_bt,
        #                     "\n"; color = :yellow)
        #     end
        # else
            trim_bt = s["base_trim"][1]
        # end
        # if s["base_trim"][2] > s["meta_btime"] + Millisecond(round(Int, t.head["npts"] * t.head["delta"] * 1e3))
        #     trim_et = s["meta_btime"] + Millisecond(round(Int, t.head["npts"] * t.head["delta"] * 1e3))
        #     if warn
        #         printstyled("End time: ", s["base_trim"][2],
        #                     " of $(s["network"]).$(s["station"]).$(s["component"]) is too late, set to: ", trim_et,
        #                     "\n"; color = :yellow)
        #     end
        # else
            trim_et = s["base_trim"][2]
        # end
        # t = trim(t, trim_bt, trim_et)
        (sbt, tw, _) = SeisTools.DataProcess.cut(t.data, SeisTools.SAC.DateTime(t.hdr), trim_bt, trim_et, 
        Millisecond(round(Int, t.hdr["delta"]*1000)))
        s["base_begintime"] = sbt
        s["base_record"] = tw
        Green.load!(s, env)
        detrendandtaper!(s["base_record"])
        detrendandtaper!(s["green_fun"])
    end
    return nothing
end

function updatestationselection!(env::Setting, status::Dict{String,Any})
    if !("stationlist" in keys(status))
        slist = String[]
        for i = 1:length(env["stations"])
            tag = join([env["stations"][i]["network"], env["stations"][i]["station"]], '.')
            push!(slist, tag)
        end
        unique!(slist)
        sort!(slist)
        status["stationlist_selected"] = Set()
        status["stationlist"] = Set(slist)
    end
    return nothing
end

returnnone(x...) = nothing

gettag(x::Dict) = x["network"] * "." * x["station"]

autoscale(x, y) =
    let
        rx = (maximum(x) - minimum(x)) / 2
        ry = (maximum(y) - minimum(y)) / 2
        cx = minimum(x) + rx
        cy = minimum(y) + ry
        (cx - rx * 1.2, cx + rx * 1.2, cy - ry * 1.2, cy + ry * 1.2)
    end

function plotstationscatter(env, status)
    tags = gettag.(env["stations"])
    lon1 = Float64[]
    lat1 = Float64[]
    tag1 = String[]
    lon2 = Float64[]
    lat2 = Float64[]
    tag2 = String[]
    for s in status["stationlist"]
        idx = findfirst(tags .== s)
        if s in status["stationlist_selected"]
            push!(lon1, env["stations"][idx]["meta_lon"])
            push!(lat1, env["stations"][idx]["meta_lat"])
            push!(tag1, s)
        else
            push!(lon2, env["stations"][idx]["meta_lon"])
            push!(lat2, env["stations"][idx]["meta_lat"])
            push!(tag2, s)
        end
    end
    if ("saveplotdata" in keys(status)) && ("saveplotdatato" in keys(status)) && status["saveplotdata"]
        matwrite(status["saveplotdatato"],
                 Dict("type" => "station", "lon1" => lon1, "lon2" => lon2,
                      "lat1" => lat1, "lat2" => lat2, "tag1" => tag1, "tag2" => tag2,
                      "elon" => env["event"]["longitude"], "elat" => env["event"]["latitude"]))
        return nothing
    end
    return nothing
end

function plotstation(env, status)
    if !("current_station" in keys(status))
        @warn "Station not specified."
        return nothing
    end
    idx = findall(gettag.(env["stations"]) .== status["current_station"])
    if isempty(idx)
        @warn "Station not specified."
        return nothing
    end
    stas = env["stations"][idx]
    bt = minimum(map(x -> x["base_begintime"], stas))
    if !("refmech" in keys(status))
        status["refmech"] = round.(rand(3) .* [360, 90, 180] - [0.0, 0.0, 90.0])
    end
    m = dc2ts(status["refmech"])
    ks = keys(status)
    flag_filter = all(x -> x in ks, ("filtwave", "filterorder", "filterband"))
    flag_filter = flag_filter ? flag_filter && status["filtwave"] : false
    flag_green = "green_fun" in keys(stas[1])
    amp = "waveamp" in keys(status) ? status["waveamp"] : 1.0
    tl = Vector{Float64}[]
    tg = Vector{Float64}[]
    ww = Vector{Float64}[]
    wg = Vector{Float64}[]
    sw = Float64[]
    sg = Float64[]
    for i = 1:3
        rbt = round(stas[i]["base_begintime"] - bt, Millisecond).value * 1e-3
        w = deepcopy(stas[i]["base_record"])
        if flag_filter
            fltr = digitalfilter(Bandpass(status["filterband"][1], status["filterband"][2];
                                          fs = 1 / stas[i]["meta_dt"]), Butterworth(status["filterorder"]))
            w = filtfilt(fltr, w)
        end
        t = (0.0:length(w)-1) .* stas[i]["meta_dt"] .+ rbt
        if stas[i]["component"] == "E"
            shift = 3.0
        elseif stas[i]["component"] == "N"
            shift = 2.0
        else
            shift = 1.0
        end
        w ./= maximum(abs, w) * 2 / amp
        w .+= shift
        push!(tl, t)
        push!(ww, w)
        if flag_green
            g = stas[i]["green_fun"] * m
            if flag_filter
                if status["filterband"][2] < 0.5 / stas[i]["green_dt"]
                    fltr = digitalfilter(Bandpass(status["filterband"][1], status["filterband"][2];
                                                  fs = 1 / stas[i]["green_dt"]), Butterworth(status["filterorder"]))
                    g = filtfilt(fltr, g)
                else
                    @warn "filterband is not correct. max frequency: $(0.5/stas[i]["green_dt"])Hz"
                end
            end
            g ./= maximum(abs, g) * 2 / amp
            g .+= shift - 0.2
            push!(tg, (0:length(g)-1) .* stas[i]["green_dt"])
            push!(wg, g)
            push!(sw, shift)
            push!(sg, shift - 0.2)
        end
    end
    maxt = maximum(vcat(tl..., tg...))
    if "xlim" in keys(status) && "xlimreset" in keys(status) && !status["xlimreset"]
        rgn = (status["xlim"][1], status["xlim"][2], 0.0, 4.0)
    else
        rgn = (0.0, maxt, 0.0, 4.0)
    end
    if ("saveplotdata" in keys(status)) && ("saveplotdatato" in keys(status)) && status["saveplotdata"]
        t = Dict{String,Any}[]
        for p in stas[1]["phases"]
            td = Dict{String,Any}()
            for k in keys(p)
                if k in ("tt", "type") || contains(k, "trim")
                    td[k] = p[k]
                end
            end
            td["at"] = round(p["at"] - bt, Millisecond).value * 1e-3
            push!(t, td)
        end
        matwrite(status["saveplotdatato"],
                 Dict("type" => "wave", "tl" => tl, "ww" => ww, "tg" => tg, "wg" => wg, "sw" => sw, "sg" => sg,
                      "bt" => Dates.format(bt, "yyyy-mm-ddTHH:MM:SS.sss"), "win" => t))
        return nothing
    end
    return nothing
end

function plotspec(env, status)
    if !("current_station" in keys(status))
        @warn "Station not specified."
        return nothing
    end
    idx = findall(gettag.(env["stations"]) .== status["current_station"])
    if isempty(idx)
        @warn "Station not specified."
        return nothing
    end
    stas = env["stations"][idx]
    fs = Vector{Float64}[]
    spec = Vector{Float64}[]
    for s in stas
        if s["component"] == "E"
            shift = 3.0
        elseif s["component"] == "N"
            shift = 2.0
        else
            shift = 1.0
        end
        w = deepcopy(s["base_record"])
        w ./= maximum(abs, w)
        w = complex(w)
        fft!(w)
        f = (1:length(w)) / length(w) / s["meta_dt"]
        idx = 1:round(Int, length(f) / 2)
        tspec = abs.(w[idx])
        tspec ./= maximum(tspec)
        tspec .+= shift
        push!(fs, f[idx])
        push!(spec, tspec)
    end
    maxf = maximum(vcat(fs...))
    if ("saveplotdata" in keys(status)) && ("saveplotdatato" in keys(status)) && status["saveplotdata"]
        tbands = Matrix{Float64}(undef, 2, 0)
        tbands = Dict{String,Any}[]
        for p in stas[1]["phases"]
            for k in keys(p)
                if contains(k, "band")
                    push!(tbands, Dict{String,Any}("phase" => p["type"], "key" => k, "band" => p[k]))
                end
            end
        end
        matwrite(status["saveplotdatato"], Dict("type" => "spec", "fs" => fs, "spec" => spec, "band" => tbands))
        return nothing
    end
    return nothing
end

function outputconfig(env, status, cmd)
    t = Dict("algorithm" => env["algorithm"], "event" => env["event"], "stations" => env["stations"])
    if length(cmd) <= 1
        tpath = normpath(env["dataroot"], "conf.jld2")
    else
        tpath = normpath(env["dataroot"], cmd[2])
        if isdir(tpath)
            tpath = normpath(tpath, "conf.jld2")
        end
    end
    jldsave(tpath; env = env, status = status)
    return nothing
end

function shortstring(x)
    t = string(x)
    return length(t) > 80 ? t[1:38] * "..." * t[end-37:end] : t
end

function parsesetconfig(env, status, cmd::Vector{String})
    if !("current_station" in keys(status)) || isempty(status["current_station"])
        @warn "station not specified"
        return Int[]
    end
    opts = ArgParser.newoptionlist()
    ArgParser.addopt!(opts, "component", "all", "c")
    (pars, args) = ArgParser.parse(cmd, opts)
    if pars["component"] == "all"
        idxs = findall(gettag.(env["stations"]) .== status["current_station"])
    else
        idxs = findall((gettag.(env["stations"]) .== status["current_station"]) .&
                       map(x -> x["component"] == pars["component"], env["stations"]))
    end
    if isempty(idxs)
        @error "Component $(pars["component"]) not exist"
    end
    return (idxs, args)
end

function printstatus(env, status, cmd)
    println("\nCurrent status:\n")
    ks = Set(keys(status))
    pop!(ks, "stationlist")
    pop!(ks, "stationlist_selected")
    ks = sort(collect(ks))
    if isempty(ks)
        ml = 8
    else
        ml = maximum(length, ks)
        ml = (floor(Int, ml / 4) + 1) * 4
        for k in ks
            println("    ", " "^(ml - length(k) - 1), k, ": ", shortstring(status[k]))
        end
    end
    println("    ", " "^(ml - 8), "station:")
    if "statusformat" in keys(status) && status["statusformat"] == "long"
        staname = String[]
        dist = Float64[]
        for s in collect(status["stationlist"])
            push!(staname, s)
            id = findfirst(gettag.(env["stations"]) .== s)
            push!(dist, env["stations"][id]["base_distance"])
        end
        for i in sortperm(dist)
            s = staname[i]
            if s in status["stationlist_selected"]
                tcolor = :green
            else
                tcolor = :red
            end
            if "current_station" in keys(status) && s == status["current_station"]
                printstyled(" "^(5 + ml), "*")
                printstyled(@sprintf("%s %.2fkm", s, dist[i]), "\n"; color = tcolor)
            else
                printstyled(" "^(6 + ml), @sprintf("%s %.2fkm", s, dist[i]), "\n"; color = tcolor)
            end
        end
    else
        for s in collect(status["stationlist"])
            if s in status["stationlist_selected"]
                tcolor = :green
            else
                tcolor = :red
            end
            if "current_station" in keys(status) && s == status["current_station"]
                print(" *")
                printstyled(s; color = tcolor)
            else
                printstyled(" ", s; color = tcolor)
            end
        end
        println("")
    end
    return nothing
end

function printevent(env, status, cmd)
    println("Event Info:\n")
    ks = keys(env["event"]) |> Set |> collect |> sort
    if !isempty(ks)
        ml = maximum(length, ks)
        ml = (floor(Int, ml / 4) + 1) * 4
        for k in ks
            println("    ", " "^(ml - length(k) - 1), k, ": ", shortstring(env["event"][k]))
        end
    end
    return nothing
end

function printalgorithm(env, status, cmd)
    println("Algorithm Status:\n")
    ks = keys(env["algorithm"]) |> Set |> collect |> sort
    if isempty(ks)
        ml = 8
    else
        ml = maximum(length, ks)
        ml = (floor(Int, ml / 4) + 1) * 4
        for k in ks
            println("    ", " "^(ml - length(k) - 1), k, ": ", shortstring(env["algorithm"][k]))
        end
    end
    return nothing
end

function printstationinfo(env, status, cmd)
    if "current_station" in keys(status) && !isempty(status["current_station"])
        # idxs = findfirst(gettag.(env["stations"]) .== status["current_station"])
        (idxs, args) = parsesetconfig(env, status, cmd)
        s = env["stations"][idxs[1]]
        println("Station: ", gettag(s))
        ks = sort(collect(keys(s)))
        ml = maximum(length, ks)
        ml = (floor(Int, ml / 4) + 1) * 4
        for k in ("station", "network", "component")
            println("    ", " "^(ml - length(k) - 1), k, ": ", shortstring(s[k]))
        end
        for k in ks
            if k in ("phases", "station", "network", "component")
                continue
            end
            println("    ", " "^(ml - length(k) - 1), k, ": ", shortstring(s[k]))
        end
        println("    phases:")
        for p in s["phases"]
            pks = keys(p) |> collect |> sort
            pml = maximum(length, pks)
            pml = (floor(Int, pml / 4) + 1) * 4
            println(" "^8, " "^(pml - 5), "type: ", p["type"])
            println(" "^8, " "^(pml - 3), "at: ", p["at"])
            println(" "^8, " "^(pml - 3), "tt: ", p["tt"])
            for pk in pks
                if pk in ("type", "at", "tt")
                    continue
                end
                println(" "^8, " "^(pml - length(pk) - 1), pk, ": ", shortstring(p[pk]))
            end
            println("")
        end
    else
        @warn "Station not specified"
    end
    return nothing
end

function setstatus!(env, status, cmd)
    if length(cmd) < 4
        help_set()
    else
        if cmd[3] == "current_station"
            status[cmd[3]] = cmd[4]
        elseif cmd[3] == "refmech"
            if length(cmd) < 6
                @warn "parameter not enough"
            else
                status["refmech"] = parse.(Float64, cmd[4:6])
            end
        elseif cmd[3] == "select"
            if cmd[4] in status["stationlist"]
                if length(cmd) < 5
                    t = !(cmd[4] in status["stationlist_selected"])
                else
                    t = parse(Bool, cmd[5])
                end
                if t
                    push!(status["stationlist_selected"], cmd[4])
                else
                    pop!(status["stationlist_selected"], cmd[4])
                end
            else
                @warn "station $(cmd[4]) not exist"
            end
        elseif cmd[3] == "filterorder"
            status["filterorder"] = parse(Int, cmd[4])
        elseif cmd[3] == "filtwave"
            status["filtwave"] = parse(Bool, lowercase(cmd[4]))
        elseif cmd[3] == "filterband"
            if length(cmd) < 5
                @warn "bandpass filter need 2 corner frequency"
            else
                status["filterband"] = parse.(Float64, cmd[4:5])
            end
        elseif cmd[3] == "xlimreset"
            status["xlimreset"] = parse(Bool, lowercase(cmd[4]))
        elseif cmd[3] == "xlim"
            if length(cmd) < 5
                @warn "bandpass filter need 2 corner frequency"
            else
                status["xlim"] = parse.(Float64, cmd[4:5])
            end
        elseif cmd[3] == "statusformat"
            status["statusformat"] = cmd[4]
        elseif cmd[3] == "waveamp"
            status["waveamp"] = parse(Float64, cmd[4])
        else
            @warn "key $key illegal."
        end
    end
    status["fresh"] = false
    printstatus(env, status, cmd)
    return nothing
end

function setalgorithm!(env, status, cmd)
    if length(cmd) < 4
        @warn "value not exist."
        return nothing
    end
    if cmd[3] == "searchdepth"
        env["algorithm"]["searchdepth"] = parse(Float64, cmd[4])
    elseif cmd[3] == "weight"
        env["algorithm"]["weight"] = parse.(Float64, cmd[4:end])
    elseif cmd[3] == "misfit"
        env["algorithm"]["misfit"] = String.(cmd[4:end])
    else
        @warn "key: $(cmd[3]) illegal"
    end
    if cmd[3] == "searchdepth"
        status["fresh"] = true
    end
    return nothing
end

function setevent!(env, status, cmd)
    if length(cmd) < 4
        @warn "value not exist."
        return nothing
    end
    if cmd[3] in ("depth", "latitude", "longitude", "magnitude")
        env["event"][cmd[3]] = parse(Float64, cmd[4])
    elseif cmd[3] == "origintime"
        env["event"]["origintime"] = DateTime(cmd[4], "y-m-dTH:M:S.s")
    else
        @warn "key: $(cmd[3]) illegal"
    end
    if cmd[3] in ("depth", "latitude", "longitude", "origintime")
        status["fresh"] = true
    end
    return nothing
end

function setstation!(env, status, cmd)
    (idxs, args) = parsesetconfig(env, status, cmd)
    if length(args) < 4
        @warn "value not exist."
        return nothing
    end
    if !("current_station" in keys(status)) || isempty(status["current_station"])
        @warn "station not specified"
        return nothing
    end
    # idxs = findall(gettag.(env["stations"]) .== status["current_station"])
    if args[3] in ("base_azimuth", "base_distance", "green_dt", "green_xl")
        t = parse(Float64, args[4])
        for i in idxs
            env["stations"][i][args[3]] = t
        end
    elseif args[3] in ("green_tsource",)
        t = parse(Float64, args[4])
        for s in env["stations"]
            s[args[3]] = t
        end
    elseif args[3] == "green_m"
        t = parse(Int, args[4])
        for i in idxs
            env["stations"][i][args[3]] = t
        end
    elseif args[3] == "base_trim"
        t = DateTime.(args[4:5])
        for i in idxs
            env["stations"][i][args[3]] = t
        end
    elseif args[3] == "green_model"
        for i in idxs
            env["stations"][i][args[3]] = args[4]
        end
    end
    status["fresh"] = true
    return nothing
end

function setphase!(env, status, cmd)
    (idxs, args) = parsesetconfig(env, status, cmd)
    if !("current_station" in keys(status)) || isempty(status["current_station"])
        @warn "station not specified"
        return nothing
    end
    # idxs = findall(gettag.(env["stations"]) .== status["current_station"])
    if length(args) == 4
        if args[3] == "type"
            if !(args[4] in map(x -> x["type"], env["stations"][idxs[1]]["phases"]))
                for i in idxs
                    push!(env["stations"][i]["phases"], Dict{String,Any}("type" => args[4]))
                end
            end
            status["fresh"] = true
            return nothing
        end
    end
    if length(args) < 5
        @warn "value not exist."
        return nothing
    end
    pidx = findfirst(x -> x["type"] == args[3], env["stations"][1]["phases"])
    if isnothing(pidx)
        @warn "phase: $(args[3]) not exist. use: set phase type typename to add new phase"
        return nothing
    end
    if args[4] in ("tt", "dtw_dt", "dtw_klim", "dtw_maxlag", "xcorr_dt", "xcorr_maxlag")
        t = parse(Float64, args[5])
        for i in idxs
            pidx = findfirst(x -> x["type"] == args[3], env["stations"][i]["phases"])
            env["stations"][i]["phases"][pidx][args[4]] = t
        end
    elseif args[4] == "at"
        t = DateTime(args[5], "y-m-dTH:M:S.s")
        for i in idxs
            pidx = findfirst(x -> x["type"] == args[3], env["stations"][i]["phases"])
            env["stations"][i]["phases"][pidx][args[4]] = t
        end
    elseif args[4] in ("dtw_order", "xcorr_order")
        t = parse(Int, args[5])
        for i in idxs
            pidx = findfirst(x -> x["type"] == args[3], env["stations"][i]["phases"])
            env["stations"][i]["phases"][pidx][args[4]] = t
        end
    elseif args[4] == "polarity_obs"
        t = parse(Float64, args[5])
        t = sign(t)
        for i in idxs
            pidx = findfirst(x -> x["type"] == args[3], env["stations"][i]["phases"])
            env["stations"][i]["phases"][pidx][args[4]] = t
        end
    elseif args[4] in ("dtw_band", "dtw_trim", "polarity_trim", "psr_trimp", "psr_trims", "xcorr_band", "xcorr_trim")
        if length(args) < 6
            @warn "need two float"
            return nothing
        else
            t = parse.(Float64, args[5:6])
            for i in idxs
                pidx = findfirst(x -> x["type"] == args[3], env["stations"][i]["phases"])
                env["stations"][i]["phases"][pidx][args[4]] = t
            end
        end
    end
    return nothing
end

function cmd_set(cmd::Vector{String})
    if length(cmd) < 2
        help_set()
        return returnnone
    end
    t = lowercase(cmd[2])
    if t == "status"
        return setstatus!
    elseif t == "algorithm"
        return setalgorithm!
    elseif t == "event"
        return setevent!
    elseif t == "station"
        return setstation!
    elseif t == "phase"
        return setphase!
    else
        @warn "key $t not exist"
    end
    return returnnone
end

function cmd_print(cmd::Vector{String})
    if length(cmd) < 2
        help_print()
        return returnnone
    end
    t = lowercase(cmd[2])
    if t == "status"
        return printstatus
    elseif t == "station"
        return printstationinfo
    elseif t == "event"
        return printevent
    elseif t == "algorithm"
        return printalgorithm
    end
    return returnnone
end

function cmd_plot(cmd::Vector{String})
    if length(cmd) < 2
        help_plot()
        return returnnone
    else
        t = lowercase(cmd[2])
        if t == "station"
            return plotstationscatter
        elseif t == "wave"
            return plotstation
        elseif t == "spec"
            return plotspec
        else
            return returnnone
        end
    end
end

function cmd_help(cmd::Vector{String})
    if length(cmd) <= 1
        help_print()
        help_plot()
        help_set()
        help_run()
        help_shell()
        help_save()
    elseif length(cmd) >= 2
        tc = lowercase(cmd[2])
        if tc == "print"
            help_print()
        elseif tc == "plot"
            help_plot()
        elseif tc == "set"
            help_set()
        elseif tc == "run"
            help_run()
        elseif tc == "save"
            help_save()
        else
            cmd_help([""])
        end
    end
    return nothing
end

function help_set()
    println("""
    Set status
            current_station String
            refmech         [Float, Float, Float]
            filtwave        Bool
            filterorder     Int
            filterband      [Float, Float]
            xlim            [Float, Float]
            xlimreset       Bool
            select stationname [Bool]
            statusformat    String
            waveamp         Float

        algorithm
        event
        station
        phase""")
    return nothing
end

function help_print()
    println("""
    prInt   status
                print current status

            station
                print infomation of current station""")
    return nothing
end

function help_plot()
    println("""
    Plot    station
                plot station scatter

            wave
                plot waveform and greenfun of current station

            spectrum
                plot frequency spectrum of current station""")
    return nothing
end

function help_run()
    println("""
    Run cmd
        run julia command""")
    return nothing
end

function help_shell()
    println("""
    shell cmd
        run command in shell env""")
    return nothing
end

function help_save()
    println("""
    save path   save configuration to specified path""")
    return nothing
end

function extractabbr(cmd::Vector{String})
    if cmd == ["p"]
        return ["plot", "wave"]
    elseif cmd[1] == "csta" && length(cmd) >= 2
        return ["set", "status", "current_station", cmd[2]]
    elseif cmd == ["status"]
        return ["print", "status"]
    elseif cmd[1] == "xlim"
        return [["set", "status"]; cmd]
    elseif lowercase(cmd[1]) == "p" && length(cmd) > 1
        return ["plot"; cmd[2:end]]
    elseif lowercase(cmd[1]) == "i"
        t = deepcopy(cmd)
        t[1] = "print"
        return t
    elseif lowercase(cmd[1]) == "s"
        t = deepcopy(cmd)
        t[1] = "set"
        return t
    elseif lowercase(cmd[1]) == "r"
        t = deepcopy(cmd)
        t[1] = "run"
        return t
    else
        return cmd
    end
end

function interpret!(inio::IO, env::Setting, status; misfits::Vector{Module} = Module[])
    if inio == stdin
        printstyled("IAT> "; color = :green)
    end
    l = String.(split(strip(readline(inio)); keepempty = false))
    if inio != stdin
        println("Command: ", join(l, ' '))
    end
    if isempty(l)
        return false
    end
    l = extractabbr(l)
    cmd = lowercase(l[1])
    if cmd in ("quit", "q", "exit")
        if inio != stdin
            print(inio, '\0')
        end
        return true
    elseif cmd == "print"
        func = cmd_print(l)
        func(env, status, l)
    elseif cmd == "plot"
        func = cmd_plot(l)
        func(env, status)
    elseif cmd == "set"
        func! = cmd_set(l)
        func!(env, status, l)
        if status["fresh"]
            fresh!(env, misfits)
        end
        status["fresh"] = false
    elseif cmd == "run"
        eval(Meta.parse(join(l[2:end], ' ')))
    elseif cmd == "shell"
        run(Cmd(l[2:end]); wait = false)
    elseif cmd == "save"
        outputconfig(env, status, l)
    elseif cmd == "help"
        cmd_help(l)
    elseif cmd == "tmpfile"
        if inio != stdin
            print(inio, abspath(status["saveplotdatato"]))
        end
    elseif cmd == "gettmp"
        if inio != stdin
            nf = filesize(status["saveplotdatato"])
            print(inio, @sprintf("%d\0", Int(nf)))
            open(status["saveplotdatato"]) do tmpio
                for _ = 1:nf
                    write(inio, read(tmpio, UInt8))
                end
            end
        else
            @warn "can't export to terminal"
            return false
        end
    else
        @warn "command: $cmd not exist"
    end
    println("<")
    if inio != stdin
        print(inio, '\0')
    end
    return false
end

function repl!(inio::IO, env::Setting, status::Dict; misfits::Vector{Module} = Module[])
    while isopen(inio)
        updatestationselection!(env, status)
        try
            t = interpret!(inio, env, status; misfits = misfits)
            if t
                break
            end
        catch err
            @warn err
            if inio != stdin
                print(inio, '\0')
            end
        end
    end
    return nothing
end

function repl!(env::Setting, status; misfits::Vector{Module} = Module[])
    repl!(stdin, env, status; misfits = misfits)
    return nothing
end

end
