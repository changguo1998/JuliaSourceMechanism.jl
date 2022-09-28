using JuliaSourceMechanism, JLD2, LinearAlgebra, Printf, CairoMakie, SeisTools

if !@isdefined result
    fname = ARGS[1]
    (mech, result) = let
        t = JLD2.load(joinpath(@__DIR__, fname))
        (t["mech"], t["result"])
    end
end

stationnames = collect(keys(result))
dists = map(x -> result[x]["dist"], stationnames)
distperm = sortperm(dists)
nstation = length(stationnames)
WIDTH = 16
HEIGHT = nstation * WIDTH / 6 * 0.3 + 2

function linearscale(v; xlim::Tuple{Float64,Float64}=(0.0, 1.0), ylim::Tuple{Float64,Float64}=(-0.5, 0.5))
    x = range(xlim[1], xlim[2]; length=length(v))
    vmin = minimum(v)
    vmax = maximum(v)
    y = @. (v - vmin) / (vmax - vmin) * (ylim[2] - ylim[1]) + ylim[1]
    return (x, y)
end

function plansubplot(subplotsize::Tuple{Int,Int}=(8, 6), xrange::Tuple{Float64,Float64}=(0.1, 0.95),
                     xsplit::Float64=0.05, ysplit::Float64=0.01, proportion::Vector{Float64}=ones(6))
    @assert subplotsize[2] == length(proportion)
    frameinfo = Matrix{Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}}(undef, subplotsize)
    useful_width = xrange[2] - xrange[1] - 5 * xsplit
    k = useful_width / sum(proportion)
    for i = 1:subplotsize[1], j = 1:subplotsize[2]
        t = xrange[1] + k * sum(proportion[1:j]) + (j - 1) * xsplit
        frameinfo[i, j] = ((t - proportion[j] * k, t),
                           (subplotsize[1] - i + 0.5 + ysplit / 2, subplotsize[1] - i + 1.5 - ysplit / 2))
    end
    return frameinfo
end

frame = plansubplot((nstation, 6), (0.12, 0.98), 0.01, 0.2, ones(6))

hideaxis = (xlabelvisible=false, xgridvisible=false, xticklabelsvisible=false, xticksvisible=false,
            ylabelvisible=false, ygridvisible=false, yticklabelsvisible=false, yticksvisible=false,
           leftspinevisible=false, rightspinevisible=false,
           topspinevisible=false, bottomspinevisible=false)

fig = Figure(; resolution=round.(Int, (WIDTH, HEIGHT) .* 100))
ax1 = Axis(fig[1, 1]; title=@sprintf("%d/%d/%d", round(Int, mech[1]), round(Int, mech[2]), round(Int, mech[3])),
           aspect=1, hideaxis...)
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[2, 1:2]; hideaxis...)

for n = 1:nstation, j = 1:6
    i = distperm[n]
    s = stationnames[i]
    p = (j > 3) ? "S" : "P"
    if mod(j - 1, 3) == 0
        c = "E"
    elseif mod(j - 1, 3) == 1
        c = "N"
    else
        c = "Z"
    end
    dt = result[s][c][p]["xcorr_dt"]
    sh = result[s][c][p]["xcorr_shift"]
    nsp = round(Int, sh / dt)
    if sh > 0
        wr = result[s][c][p]["xcorr_rec"][nsp:end]
        ws = result[s][c][p]["xcorr_syn"][1:end-nsp+1]
    else
        wr = result[s][c][p]["xcorr_rec"][1:end+nsp]
        ws = result[s][c][p]["xcorr_syn"][1-nsp:end]
    end
    (tx, tr) = linearscale(wr; xlim=frame[n, j][1], ylim=frame[n, j][2])
    lines!(ax3, tx, tr; color=:black)
    (tx, tr) = linearscale(ws; xlim=frame[n, j][1], ylim=frame[n, j][2])
    lines!(ax3, tx, tr; color=:red)
    hl = (frame[n, j][2][1] + frame[n, j][2][2]) / 2
    text!(ax3, 0.12, hl; text=stationnames[i], align=(:right, :baseline))
end

heatmap!(ax1, range(; start=-1.0, stop=1.0, length=1001), range(; start=-1.0, stop=1.0, length=1001),
         SeisTools.Source.beachball_bitmap(SeisTools.Source.MomentTensor(mech[1], mech[2], mech[3]);
                                           resolution=(1001, 1001)) |>
         permutedims .|> sign;
         colormap=:binary)
l = SeisTools.Source.beachball_sdrline(SeisTools.Source.MomentTensor(mech[1], mech[2], mech[3]))
lines!(ax1, map(x->(x[2], x[1]), l.l1), color=:black, width=2)
lines!(ax1, map(x->(x[2], x[1]), l.l2), color=:black, width=2)
lines!(ax1, map(x->(x[2], x[1]), l.edge), color=:black, width=2)
rowsize!(fig.layout, 1, Fixed(round(Int, WIDTH*0.25*0.618*100)))
save("test.png", fig, px_per_unit=4)
