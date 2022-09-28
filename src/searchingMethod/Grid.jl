module Grid

STEP = [5.0, 5.0, 5.0]

function setstep!(strike::Real, dip::Real, rake::Real)
    global STEP
    STEP[1] = strike
    STEP[2] = dip
    STEP[3] = rake
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
    for s = 0.0:STEP[1]:355.0, d = 0.0:STEP[2]:90.0, r = -90.0:STEP[3]:90.0
        push!(newsdr, [s, d, r])
    end
    return newsdr
end
end
