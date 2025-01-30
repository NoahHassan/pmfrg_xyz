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

function ChiHom(T, w, j)
    b = 1 / T
    return (
        -(exp(b) - 1 - b) / (2*(exp(b) + 3)) * (1 / (1 + w^2))
    )
end

Trange = 0.5:0.25:3.0
x = range(0.5, 3.0, length=100)

chiVals = []

using JLD2
for T in Trange
    sol = load_object("Tflow/flow$T.jld2")
    chi = sol[end][2]
    append!(chiVals, chi)
end

fig = Figure()
ax = Axis(fig[1,1], ylabel = L"χ",xlabel = L"T", title = "Dimer Jx = 1.0, Jy = -0.6, Jz = 0.7")
lines!(ax, x, ChiExact12.(x, 0, 1.0, -0.6, 0.7), color=:black, linewidth=2)
scatter!(ax,Trange,chiVals, marker=:diamond, color=:red, markersize=10)
# lines!(ax, x, ChiHom.(x, 0, 1.0))
# axislegend(ax, position=:rt)
display("image/png", fig)