using ArgumentProcessor, JuliaSourceMechanism, TOML, Dates, DelimitedFiles, Printf, JLD2, FFTW, Statistics, DSP, Sockets

# a user defined initial function
include("lib.jl")

# set misfit model
misfit_module = [XCorr, Polarity]

addopt!("port"; abbr = "P", fmt = " %d", default=" 12345", help = "Port number opened to GUI program")
addopt!("dataroot"; abbr = "R", fmt=" %s", default = " "*pwd(), help = "Root directory of event")
addopt!("jldfile"; abbr="J", fmt=" %s", default=" ", help="Data file")
INPUTS = ArgumentProcessor.parse(ARGS)

if isdir(INPUTS.dataroot)
    dataroot = abspath(INPUTS.dataroot)
else
    if isdir(normpath(pwd(), dataroot))
        dataroot = abspath(pwd(), dataroot)
    else
        @error "dataroot not exist"
    end
end
if isfile(normpath(INPUTS.dataroot, INPUTS.jldfile))
    (env, status) = let
        t = load(normpath(INPUTS.dataroot, INPUTS.jldfile))
        (t["env"], t["status"])
    end
else
    env = loadandinit(misfit_module, dataroot)
    status = Dict{String, Any}()
end
status["saveplotdata"] = true
status["saveplotdatato"] = abspath(".tmpplot.mat")

@info "Event dir: $dataroot\nTmp plot file: $(status["saveplotdatato"])"

# uncomment codes below while using matlab interface
svr = listen(IPv4(0), INPUTS.port)
@info "Listening localhost:$(INPUTS.port)\nwaiting for connection..."
sock = accept(svr)
@info "Connection established"

JuliaSourceMechanism.InteractiveTest.repl!(sock, env, status; misfits = misfit_module)

# uncomment codes below while using CLI
# JuliaSourceMechanism.InteractiveTest.repl!(env, status; misfits = misfit_module)
