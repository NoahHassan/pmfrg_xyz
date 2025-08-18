import Pkg
Pkg.activate(".")

using CairoMakie

function ChiExact12(T, w, jx, jy, jz)
    b = 1 / T
    Jp = (jy + jx) / 2.0
    Jm = (jx - jy) / 2.0
    Jyz = (jz + jy) / 2.0
    Jxz = (jz + jx) / 2.0

    v = 0
    if(Jm == 0)
        v = exp(-10)
    end

    Z = 2.0 * (1.0 + exp(-b * Jp) + exp(-b * Jxz) + exp(-b * Jyz))
    A = -(1.0 - exp(-b * Jp)) / Jp
    B = -(1.0 - exp(b * Jm)) * exp(-b * Jxz) / Jm

    if(isapprox(Jm, 0.0, atol=1e-5))
        B = (b + b^2 * Jm) * exp(-b * Jxz)
    end

    result = (A + B) / Z

    return result
end

function ChiExact11(T, w, jx, jy, jz)
    b = 1 / T
    Jp = (jy + jx) / 2.0
    Jm = (jx - jy) / 2.0
    Jyz = (jz + jy) / 2.0
    Jxz = (jz + jx) / 2.0

    v = 0
    if(Jm == 0)
        v = exp(-10)
    end

    Z = 2.0 * (1.0 + exp(-b * Jp) + exp(-b * Jxz) + exp(-b * Jyz))
    A = (1.0 - exp(-b * Jp)) / Jp
    B = -(1.0 - exp(b * Jm)) * exp(-b * Jxz) / Jm

    if(isapprox(Jm, 0.0, atol=1e-5))
        B = (b + b^2 * Jm) * exp(-b * Jxz)
    end

    result = (A + B) / Z

    return result
end

Chi_z(T, jx, jy, jz) = ChiExact12(T, 0.0, jx, jy, jz)
Chi_x(T, jx, jy, jz) = ChiExact12(T, 0.0, jz, jy, jx)
Chi_y(T, jx, jy, jz) = ChiExact12(T, 0.0, jx, jz, jy)
ChiLocal_z(T, jx, jy, jz) = ChiExact11(T, 0.0, jx, jy, jz)
ChiLocal_x(T, jx, jy, jz) = ChiExact11(T, 0.0, jz, jy, jx)
ChiLocal_y(T, jx, jy, jz) = ChiExact11(T, 0.0, jx, jz, jy)

function ChiHom(T, w, j)
    b = 1 / T
    return (
        -(exp(b) - 1 - b) / (2*(exp(b) + 3)) * (1 / (1 + w^2))
    )
end

chiVals_x = []
chiVals_y = []
chiVals_z = []

using JLD2
data = load_object("ChiDimer.jld2")

length(data)

for n in eachindex(data)
    append!(chiVals_x, data[n][1].Chi_x[2])
    append!(chiVals_y, data[n][1].Chi_y[2])
    append!(chiVals_z, data[n][1].Chi_z[2])
end

isotropy = [0.4, 0.0, 1.0]

T = [data[n][end] for n in eachindex(data)]

fig = Figure()
ax = Axis(
    fig[1,1],
    ylabel = L"χ",
    xlabel = L"T",
    title = "Jx = 0.4, Jy = 0.0, Jz = 1.0"
    )
n = 350
lines!(ax, T[n:end], Chi_x.(T[n:end], isotropy[1], isotropy[2], isotropy[3]), linewidth=2, label=L"χ_x")
lines!(ax, T[n:end], Chi_y.(T[n:end], isotropy[1], isotropy[2], isotropy[3]), linewidth=2, label=L"χ_y")
lines!(ax, T[n:end], Chi_z.(T[n:end], isotropy[1], isotropy[2], isotropy[3]), linewidth=2, label=L"χ_z")
scatter!(ax, T[n:end], chiVals_x[n:end], marker=:diamond, markersize=5)
scatter!(ax, T[n:end], chiVals_y[n:end], marker=:diamond, markersize=5)
scatter!(ax, T[n:end], chiVals_z[n:end], marker=:diamond, markersize=5)
# lines!(ax, x, ChiHom.(x, 0, 1.0))
axislegend(ax, position=:rt)
display("image/png", fig)



### Check self energy error

using JLD2;

function iSigmaExact(T, w, jx, jy, jz)
    b = 1 / T
    Ti = T * pi
    Jp = (jy + jx) / 2.0
    Jm = (jx - jy) / 2.0
    Jyz = (jz + jy) / 2.0
    Jxz = (jz + jx) / 2.0

    v = 0
    if(Jm == 0)
        v = exp(-10)
    end

    Z = (1.0 + exp(-b * Jp) + exp(-b * Jxz) + exp(-b * Jyz))
    A = (1.0 + exp(-b * Jp)) * Ti / (Jp^2 + Ti^2)
    B = (1.0 + exp(b * Jm)) * exp(-b * Jxz) * Ti / (Jm^2 + Ti^2)

    cornelator = -(A + B) / Z

    return -(Ti + 1.0 / cornelator)
end

iSig_z(T, jx, jy, jz) = iSigmaExact(T, 0.0, jx, jy, jz)
iSig_x(T, jx, jy, jz) = iSigmaExact(T, 0.0, jz, jy, jx)
iSig_y(T, jx, jy, jz) = iSigmaExact(T, 0.0, jx, jz, jy)

Trange = exp10.(range(-1, 1, length=20))
x = range(0.1, 10.00, length=100)
isotropy = [-1.0, -1.0, 1.0]

sig_z = Float64[]
sig_x = Float64[]
sig_y = Float64[]

for i in axes(Trange,1)
    sol = load_object("TflowSigma4/flow$i.jld2")
    append!(sig_z, abs(iSig_z(Trange[i], isotropy[1], isotropy[2], isotropy[3]) - sol[end].x[4][1]))
    append!(sig_x, abs(iSig_x(Trange[i], isotropy[1], isotropy[2], isotropy[3]) - sol[end].x[2][1]))
    append!(sig_y, abs(iSig_y(Trange[i], isotropy[1], isotropy[2], isotropy[3]) - sol[end].x[3][1]))
end

function opt_a(xVals, yVals)
    a = 0.0
    b = 0.0
    for i in axes(xVals, 1)
        a += yVals[i] / xVals[i]^3
        b += 1.0 / xVals[i]^6
    end

    return a / b
end

function errFunc(a, x)
    return a / x^3
end

a_z = opt_a(collect(Trange), sig_z)
a_y = opt_a(collect(Trange), sig_y)
a_x = opt_a(collect(Trange), sig_x)

fig = Figure()
ax = Axis(fig[1,1], ylabel = L"iΣ^y_{ex} - iΣ^y",xlabel = L"T", title = "Jx = -1.0, Jy = -1.0, Jz = 1.0", xscale=log, yscale=log)
# lines!(ax, x, errFunc.(a_x, x), linewidth=2, label=L"iΣ_x")
# lines!(ax, x, errFunc.(a_y, x), linewidth=2, label=L"iΣ_y")
lines!(ax, x, errFunc.(a_y, x), linewidth=2, label=L"1/T^3", color=:black)
# scatter!(ax, Trange, sig_x, marker=:diamond, markersize=10, label="x")
# scatter!(ax, Trange, sig_y, marker=:diamond, markersize=10, label="y")
scatter!(ax, Trange, sig_y, marker=:diamond, markersize=10, label=L"iΣ^y_{ex} - iΣ^y", color=:red)
axislegend(ax, position=:rt)
display("image/png", fig)

### T flow check ###
chi, tri = load_object("dimer_chi.jld2")
### StateType: f_int, iSigma (3 elements), Gamma

isotropy = [0.1, -0.6, 1.0]
fig = Figure()
ax = Axis(fig[1,1], title="T-flow [0.1, -0.6, 1.0]")
scatter!(ax, exp.(tri), getindex.(chi,1))
scatter!(ax, exp.(tri), getindex.(chi,2))
lines!(ax, exp.(tri), ChiLocal_z.(exp.(tri), isotropy[1], isotropy[2], isotropy[3]), linewidth=2, label="local")
lines!(ax, exp.(tri), Chi_z.(exp.(tri), isotropy[1], isotropy[2], isotropy[3]), linewidth=2, label="on site")
axislegend(ax, position=:rt)
display("image/png", fig)