module VelocityModel

export readmodel, readmodel_crust10

function _get_ilat_ilon(lat::Real, lon::Real)
    ilat = floor(Int, 90.5-lat) + 1
    ilat -= div(ilat, 180)
    ilon = floor(Int, lon+180.0) + 1
    ilon -= div(ilon, 360)
    return (ilat, ilon)
end

function _read_crust_1(ilat::Int, ilon::Int)
    io = open(joinpath(@__DIR__, "dat", "crust1.bin"))
    seek(io, ((ilat-1)*360+ilon-1)*144)
    v = zeros(Float32, 9, 4)
    read!(io, v)
    return v
end

function _recal_model(m::Matrix{Float32})
    thick = -diff(m[:, 1])
    idxs = findall(>(0.0), thick)
    push!(idxs, 9)
    return m[idxs, :]
end

function readmodel_crust10(lat::Real, lon::Real)
    (ilat, ilon) = _get_ilat_ilon(lat, lon)
    m = _recal_model(_read_crust_1(ilat, ilon))
    m[:, 1] .*= -1
    return round.(m, digits=2)
end

function readmodel(lat::Real, lon::Real, model::Union{String,Symbol}=:crust10)
    return readmodel_crust10(lat, lon)
end

end