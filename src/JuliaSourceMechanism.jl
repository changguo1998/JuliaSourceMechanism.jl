module JuliaSourceMechanism
using Dates, DelimitedFiles, DSP, FFTW, LinearAlgebra, Printf, Statistics, TOML, SeisTools

Setting = Dict{String,Any}

export Green, XCorr, Polarity, PSR, DTW, AbsShift, RelShift, CAP,
       Grid,
       Setting, dc2ts, sacDateTime, trim, preprocess!, inverse!, CAPmethod!, parsesac,
       mergemeta, checkconfiguration, buildstationconfiguration

# * basic macros and functions

"""
    @must(cond, text = "")

like `@assert`, but is user defined and will always be executed.
if `cond` is false, the macro will throw an error with `text`
"""
macro must(cond, text = "")
    return :(if !($(esc(cond)))
                 error($(esc(text)))
             end)
end

"""
    @hadbetter(cond, text = "")

like `@must`, but only warning when `cond` is false with information `text`
"""
macro hadbetter(cond, text = "")
    return :(if !($(esc(cond)))
                 @warn($text)
             end)
end

# * modules
# include("Seis.jl")
include("mathematics.jl")
include("system.jl")
include("Green.jl")
for m in filter(endswith(".jl"), readdir(normpath(@__DIR__, "misfits/"); join = true))
    include(m)
end
for s in filter(endswith(".jl"), readdir(normpath(@__DIR__, "searchingMethod/"); join = true))
    include(s)
end
# include("InteractiveTest.jl")
end # module
