module Grid

START = [0.0, 0.0, -90.0]
STEP = [5.0, 5.0, 5.0]
STOP = [355.0, 90.0, 90.0]

function setstart!(strike::Real, dip::Real, rake::Real)
    global START
    START[1] = strike
    START[2] = dip
    START[3] = rake
end

function setstep!(strike::Real, dip::Real, rake::Real)
    global STEP
    STEP[1] = strike
    STEP[2] = dip
    STEP[3] = rake
end

function setstop!(strike::Real, dip::Real, rake::Real)
    global STOP
    STOP[1] = strike
    STOP[2] = dip
    STOP[3] = rake
end

function continueloop!(sdr::Vector{Vector{Float64}}, misfit::Vector{Float64}, status::Dict, env::Dict)
    if isempty(sdr) && isempty(misfit) && !("grid_itr" ∈ keys(status))
        return true
    else
        return false
    end
end

function newparameters(sdr::Vector{Vector{Float64}}, misfit::Vector{Float64})
    global STEP
    newsdr = Vector{Float64}[]
    for s = START[1]:STEP[1]:STOP[1], d = START[2]:STEP[2]:STOP[2], r = START[3]:STEP[3]:STOP[3]
        push!(newsdr, [s, d, r])
    end
    return newsdr
end
end
