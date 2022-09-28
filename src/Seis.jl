module Seis
using Statistics, Printf

abstract type Frame <: Any end

struct WaveFrame <: Frame
    type::AbstractString
    head::Dict
    data::AbstractArray
    function WaveFrame(type = "default", head = Dict(), data = [])
        return new(type, head, data)
    end
end

struct SACFrame <: Frame
    head::Dict
    data::Array{Array{Float32,1},1}
    function SACFrame(head, data)
        return new(head, data)
    end
end

struct SEGYFrame <: Frame
    head::Dict
    data::Vector
    function SEGYFrame(head, data)
        return new(head, data)
    end
end

module HEADER
module tmp
sacheadlist = ("delta", "depmin", "depmax", "scale", "odelta", "b", "e", "o", "a", "internal1", "t0", "t1", "t2", "t3",
               "t4", "t5", "t6", "t7", "t8", "t9", "f", "resp0", "resp1", "resp2", "resp3", "resp4", "resp5", "resp6",
               "resp7", "resp8", "resp9", "stla", "stlo", "stel", "stdp", "evla", "evlo", "evel", "evdp", "mag",
               "user0", "user1", "user2", "user3", "user4", "user5", "user6", "user7", "user8", "user9", "dist", "az",
               "baz", "gcarc", "internal2", "internal3", "depmen", "cmpaz", "cmpinc", "xminimun", "xmaximum",
               "yminimum", "ymaximum", "unused1", "unused2", "unused3", "unused4", "unused5", "unused6", "unused7",
               "nzyear", "nzjday", "nzhour", "nzmin", "nzsec", "nzmsec", "nvhdr", "norid", "nevid", "npts", "internal4",
               "nwfid", "nxsize", "nysize", "unused8", "iftype", "idep", "iztype", "unused9", "iinst", "istreg",
               "ievreg", "ievtyp", "iqual", "isynth", "imagtyp", "imagsrc", "unused10", "unused11", "unused12",
               "unused13", "unused14", "unused15", "unused16", "unused17", "leven", "lpspol", "lovrok", "lcalda",
               "unused18", "kstnm", "kevnm", "khole", "ko", "ka", "kt0", "kt1", "kt2", "kt3", "kt4", "kt5", "kt6",
               "kt7", "kt8", "kt9", "kf", "kuser0", "kuser1", "kuser2", "kcmpnm", "knetwk", "kdatrd", "kinst")

enumeratevars = ("iftype", "idep", "iztype", "ievtyp", "iqual", "isynth", "imagtyp", "imagsrc")
logicalvars = ("leven", "lpspol", "lovrok", "lcalda", "unused18")
iftype = Dict([(1, "ITIME"), (2, "IRLIM"), (3, "IAMPH"), (4, "IXY"), (51, "IXYZ")])
idep = Dict([(5, "IUMKN"), (6, "IDISP"), (7, "IVEL"), (8, "IACC"), (50, "IVOLTS")])
iztype = Dict([(5, "IUNKN"), (9, "IB"), (10, "IDAY"), (11, "IO"), (12, "IA"), (13, "IT0"), (14, "IT1"), (15, "IT2"),
               (16, "IT3"), (17, "IT4"), (18, "IT5"), (19, "IT6"), (20, "IT7"), (21, "IT8"), (22, "IT9")])
ievtyp = Dict([(5, "IUNKN"), (37, "INUCL"), (38, "IPREN"), (39, "IPOSTN"), (40, "IQUAKE"), (41, "IPREQ"),
               (42, "IPOSTQ"), (43, "ICHEM"), (44, "IOTHER"), (70, "IQB"), (71, "IQB1"), (72, "IQB2"), (73, "IQBX"),
               (74, "IQMT"), (75, "IEQ"), (76, "IEQ1"), (77, "IEQ2"), (78, "IME"), (79, "IEX"), (80, "INU"),
               (81, "INC"), (82, "IO_"), (83, "IL"), (84, "IR"), (85, "IT"), (86, "IU")])
iqual = Dict([(44, "IOTHER"), (45, "IGOOD"), (46, "IGLCH"), (47, "IDROP"), (48, "ILOSN")])
isynth = Dict([(49, "IRLDA")])
imagtyp = Dict([(52, "IMB"), (53, "IMS"), (54, "IML"), (55, "IMW"), (56, "IMD"), (57, "IMX")])
imagsrc = Dict([(58, "INEIC"), (59, "IPDE"), (60, "IISC"), (61, "IREB"), (62, "IUSGS"), (63, "IBRK"), (64, "ICALTECH"),
                (65, "ILLNL"), (66, "IEVLOC"), (67, "IJSOP"), (68, "IUSER"), (69, "IUNKNOWN")])

iftyper = Dict([(i.second, i.first) for i in iftype])
idepr = Dict([(i.second, i.first) for i in idep])
iztyper = Dict([(i.second, i.first) for i in iztype])
ievtypr = Dict([(i.second, i.first) for i in ievtyp])
iqualr = Dict([(i.second, i.first) for i in iqual])
isynthr = Dict([(i.second, i.first) for i in isynth])
imagtypr = Dict([(i.second, i.first) for i in imagtyp])
imagsrcr = Dict([(i.second, i.first) for i in imagsrc])
sacheadtrans = (int2other = Dict([("iftype", iftype), ("idep", idep), ("iztype", iztype), ("ievtyp", ievtyp),
                                  ("iqual", iqual), ("isynth", isynth), ("imagtyp", imagtyp), ("imagsrc", imagsrc)]),
                other2int = Dict([("iftype", iftyper), ("idep", idepr), ("iztype", iztyper), ("ievtyp", ievtypr),
                                  ("iqual", iqualr), ("isynth", isynthr), ("imagtyp", imagtypr), ("imagsrc", imagsrcr)]))
segyheadlist = ()
end
sac = (headlist = tmp.sacheadlist, enumeratevars = tmp.enumeratevars, headtrans = tmp.sacheadtrans,
       logicalvars = tmp.logicalvars)
segy = ()
end

export SACFrame
# sac related functions

function simplifystring(x::Vector{T}) where {T<:Union{Char,UInt8}}
    t = Char.(x)
    idx = findfirst(t .== '\0')
    if isnothing(idx)
        return String(strip(String(t), [' ', '\0']))
    else
        return String(strip(String(t[1:idx-1]), [' ', '\0']))
    end
end

"""
# readsachead( path::IO )

read head of sac binary file
"""
function readsachead(io::IO)
    hf = zeros(Float32, 70)
    hi = zeros(Int32, 40)
    hc = zeros(UInt8, 192)
    read!(io, hf)
    read!(io, hi)
    read!(io, hc)
    # change undefined
    htu               = []
    hf                = Float64.(hf)
    hf[hf.==-12345.0] .= NaN
    hi                = Float64.(hi)
    hi[hi.==-12345.0] .= NaN
    varname           = HEADER.sac.headlist
    # construct head
    for i = 1:70
        htu = [htu; (varname[i], hf[i])]
    end
    for i = 1:40
        htu = [htu; (varname[70+i], hi[i])]
    end
    htu = [htu; (varname[111], simplifystring(hc[1:8]))]
    htu = [htu; (varname[112], simplifystring(hc[9:24]))]
    for i = 3:23
        htu = [htu
               (varname[110+i], simplifystring(hc[(i*8+1):(i*8+8)]))]
    end
    head = Dict(htu)
    # numerator and logical
    for i in HEADER.sac.enumeratevars
        if !isnan(head[i])
            head[i] = HEADER.sac.headtrans.int2other[i][head[i]]
        end
    end
    for i in HEADER.sac.logicalvars
        if !isnan(head[i])
            head[i] = Bool(head[i])
        end
    end
    return head
end

function readsachead(path::AbstractString)
    return open(readsachead, path, "r")
end

"""
readsac

    read sac binary file

    =========================

    readsac(path::AbstractString)

    =========================

    readsac(io::IO)
"""
function readsac(io::IO)
    # read head
    head = readsachead(io)
    # read data
    td = zeros(Float32, Int(head["npts"]))
    if head["npts"] >= 1 && (head["iftype"] in ["IXY", "IRLIM", "IAMPH", "IXYZ", "ITIME"]) && !isnan(head["leven"])
        read!(io, td)
        data = [td]
        if !head["leven"] || (head["iftype"] in ["IAMPH", "IRLIM"])
            read!(io, td)
            data = [data, td]
        end
    else
        data = []
    end
    return SACFrame(head, data)
end

function readsac(path::AbstractString)
    return open(readsac, path, "r")
end

"""
writesac

    this function won't check if the head is right!!!

    =======================

    writesac(io::IO, data::SACFrame)

    =======================

    writesac(path:AbstractString="./seismogram.sac", data::SACFrame)
"""
# todo   add function to check header and data
# !  this function don't check right now!!
function writesac(io::IO, data::SACFrame = SACFrame(Dict(), []))
    for i = 1:110
        hname = HEADER.sac.headlist[i]
        if isequal(data.head[hname], NaN)
            if i <= 70
                write(io, Float32(-12345.0))
            else
                write(io, Int32(-12345))
            end
        else
            if i <= 70
                write(io, Float32(data.head[hname]))
            else
                if hname in HEADER.sac.enumeratevars
                    write(io, Int32(HEADER.sac.headtrans.other2int[hname][data.head[hname]]))
                else
                    write(io, Int32(data.head[hname]))
                end
            end
        end
    end
    hname = HEADER.sac.headlist[111]
    v = data.head[hname]
    for i = 1:8
        if i <= length(v)
            write(io, UInt8(v[i]))
        else
            write(io, UInt8(' '))
        end
    end
    hname = HEADER.sac.headlist[112]
    v = data.head[hname]
    for i = 1:16
        if i <= length(v)
            write(io, UInt8(v[i]))
        else
            write(io, UInt8(' '))
        end
    end
    for i = 113:length(HEADER.sac.headlist)
        hname = HEADER.sac.headlist[i]
        v = data.head[hname]
        for j = 1:8
            if j <= length(v)
                write(io, UInt8(v[j]))
            else
                write(io, UInt8(' '))
            end
        end
    end
    for i in data.data
        write(io, i)
    end
end

function writesac(path::AbstractString, data::SACFrame = SACFrame(Dict(), []))
    return open(x -> writesac(x, data), path, "w")
end

function writesac(data::SACFrame = SACFrame(Dict(), []))
    return open(x -> writesac(x, data), "./seismogram.SAC", "w")
end

"""
# newsachead()

return new sac file head
"""
function newsachead(; b::AbstractFloat = 0.0, leven::Bool = true, idep::AbstractString = "IVEL",
                    iftype::AbstractString = "ITIME", lovrok::Bool = true, nvhdr::Int = 6, lcalda::Bool = true)
    htu = []
    varname = HEADER.sac.headlist
    for i = 1:105
        htu = [htu; (varname[i], NaN32)]
    end
    for i = 106:110
        htu = [htu; (varname[i], false)]
    end
    for i = 1:23
        htu = [htu; (varname[110+i], "-12345")]
    end
    head = Dict(htu)
    head["nvhdr"] = nvhdr
    head["b"] = b
    head["leven"] = leven
    head["idep"] = idep
    head["iftype"] = iftype
    head["lovrok"] = lovrok
    head["lcalda"] = lcalda
    return head
end

function WaveFrame(s::SACFrame)
    return WaveFrame("sac", s.head, s.data)
end

function SACFrame(s::WaveFrame)
    headkeys = keys(s.head)
    h = newsachead()
    for k in HEADER.sac.headlist
        if k in headkeys
            h[k] = s.head[k]
        end
    end
    return SACFrame(h, s.data)
end

function stdname(s::SACFrame, type::Int = 1, qtag::AbstractChar = 'D')
    if type == 1
        return @sprintf("%04d.%03d.%02d.%02d.%02d.%04d.%s.%s.%s.%s.%s.SAC", s.head["nzyear"], s.head["nzjday"],
                        s.head["nzhour"], s.head["nzmin"], s.head["nzsec"], s.head["nzmsec"], s.head["knetwk"],
                        s.head["kstnm"], s.head["khole"], s.head["kcmpnm"], qtag)
    elseif type == 2
        return @sprintf("%s.%s.%s.%s.%s.%04d.%03d.%02d%02d%02d.SAC", s.head["knetwk"], s.head["kstnm"], s.head["khole"],
                        s.head["kcmpnm"], qtag, s.head["nzyear"], s.head["nzjday"], s.head["nzhour"], s.head["nzmin"],
                        s.head["nzsec"])
    end
end

export readsac, readsachead, writesac, newsachead
end
