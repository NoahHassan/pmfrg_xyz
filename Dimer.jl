using Pkg
Pkg.activate(".")
include("pmfrg_xyz.jl")

using .pmfrg_xyz

using JLD2
using RecursiveArrayTools
using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,RecursiveArrayTools,StructArrays
using SpinFRGLattices.StaticArrays
using CairoMakie

System = getPolymer(2)
isotropy = [1.0, 0.6, 0.3]

Par = Params(
    System,
    T = 0.5,
    N = 8,
    accuracy = 1e-5,
    lambda_max = 100.,
    lambda_min = .01
)

@time sol = SolveFRG(Par, isotropy, method = DP5());
tri = LinRange(3,-2,20)
save_object("dimer_flow_noah.jld2", [(sol(t), exp(t), Par) for t in tri])

sol = load_object("dimer_flow_noah.jld2")

chiR = [getChi_z(s...) for s in sol]
fig = Figure()
ax = Axis(fig[1,1], ylabel = L"χ",xlabel = L"Λ")

scatterlines!(ax,exp.(tri),getindex.(chiR,1))
scatterlines!(ax,exp.(tri),getindex.(chiR,2))
display("image/png", fig)

## Dimer against T

System = getPolymer(2)
isotropy = [-1.0, -1.0, 1.0]

trihi = LinRange(3,-2,20)
Trange = exp10.(range(-1, 1, length=20))
print(Trange)

for n in axes(Trange,1)
    println("Solving $n...\n")
    sleep(2.0)

    Par = Params(
        System,
        T = Trange[n],
        N = 8,
        accuracy = 1e-5,
        lambda_max = 100.,
        lambda_min = .01,
        lenIntw_acc = 24
    )

    @time sol = SolveFRG(Par, isotropy, method = DP5());
    sig_z = [sol(t) for t in trihi]
    save_object("TflowSigma4/flow$n.jld2", sig_z)
end