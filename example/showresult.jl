using JLD2, Printf, DSP, LinearAlgebra

(result, mech) = let
    t = load("result.jld2")
    (t["result"], t["mech"])
end

symsign(x::Real) = let
    if x > 0.1
        return '+'
    elseif x < -0.1
        return '-'
    else
        return 'x'
    end
end

function similarity(x, y, s, dt)
    si = round(Int, s/dt)
    ax = norm(x, 2)
    ay = norm(y, 2)
    c = 0.0
    if si >= 0
        for i = 1:length(x)-si
            c += x[i+si]*y[i]
        end
    else
        for i = 1:length(x)+si
            c += x[i]*y[i-si]
        end
    end
    return c/ax/ay
end

infoline(c) = @sprintf("%c       %c    %5.1f  %5.2f  %5.1f  %5.2f", symsign(c["P"]["polarity_rec"]),
        symsign(c["P"]["polarity_syn"]), 100*similarity(c["P"]["xcorr_rec"], c["P"]["xcorr_syn"], c["P"]["xcorr_shift"], c["P"]["xcorr_dt"]),
        c["P"]["xcorr_shift"], 100*similarity(c["S"]["xcorr_rec"], c["S"]["xcorr_syn"], c["S"]["xcorr_shift"], c["S"]["xcorr_dt"]),
        c["S"]["xcorr_shift"])

println("Result:")
@printf("strike: %d  dip: %d  rake: %d\n", round.(Int, mech)...)
println("")
stations = collect(keys(result))
ds = map(t->result[t]["dist"], stations)
idx = sortperm(ds)
begin
    @printf(" station  cmp   distance  pol_rec pol_syn Pxcorr Pshift Sxcorr Sshift\n")
    for i = 1:length(idx), cmp in ("E", "N", "Z")
        s = stations[idx[i]]
        print("  ", s, "   ")
        print(cmp, "    ")
        @printf("%6.2f       ", result[s]["dist"])
        println(infoline(result[s][cmp]))
    end
end
