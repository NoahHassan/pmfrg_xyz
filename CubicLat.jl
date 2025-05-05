using Pkg
Pkg.activate(".")
include("pmfrg_xyz.jl")

using .pmfrg_xyz

using JLD2
using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,RecursiveArrayTools,StructArrays
using SpinFRGLattices.StaticArrays
using SpinFRGLattices.SimpleCubic

NLen = 2
J1 = 1
J2 = 0.0
couplings = [J1,J2]
isotropy = [1.0, 1.0, 1.0]

System = getCubic(NLen,couplings)

Par = Params(
    System,
    T=0.5,
    N = 2,
    accuracy = 1e-3,
    lambda_max = exp(10.),
    lambda_min = exp(-10.),
)

@time sol = SolveFRG(Par, isotropy, method = DP5())
save_object("LamFlow$(NLen).jld2", sol[end])

println(exp(10.))
println(exp(-10.))
## Evaluation Square lattice
@time begin
    
    using PMFRGEvaluation
    using CairoMakie #for plotting. You can use whatever plotting package you like of course

    System = SimpleCubic.getCubic(NLen)
    Lattice = LatticeInfo(System,SimpleCubic)
    let 
        L = -10
        chi_R = getChi_z(sol[end], -10, Par)
        
        chi = getFourier(chi_R, Lattice)
        
        k = LinRange(-2pi,2pi,300)

        sq2 = 1.0/sqrt(2)
        sq6 = 1.0/sqrt(6)
        
        chik = [chi(sq6 * 2.0 * x, sq6 * x + sq2 * y, sq6 * x - sq2 * y) for x in k, y in k]
        
        fig, ax, hm = heatmap(k,k,chik,axis = (;aspect = 1))
        ax.title = "Cubic lattice Λ-flow"
        Colorbar(fig[1,2],hm)
        display("image/png", fig)
    end

end

### T-flow

using Pkg
Pkg.activate(".")
include("Tpmfrg_xyz.jl")

using .Tpmfrg_xyz

using JLD2
using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,RecursiveArrayTools,StructArrays
using SpinFRGLattices.StaticArrays
using SpinFRGLattices.SimpleCubic

NLen = 4
J1 = 1.0
J2 = 0.0
System = getCubic(NLen, [J1, J2])
isotropy = [1.0, 1.0, 1.0]

Par = Params(
    System,
    N = 4,
    accuracy = 1e-3,
    temp_max = exp(10.0),
    temp_min = exp(0.1)
)

@time sol = SolveFRG(Par, isotropy, method = DP5());

using PMFRGEvaluation
using CairoMakie

System = getCubic(NLen)
Lattice = LatticeInfo(System, SimpleCubic)

let 
    T = 10.0
    chi_R = getChi_z(sol(T), T, Par)
    
    chi = getFourier(chi_R, Lattice)
    
    k = LinRange(-2pi,2pi,300)

    sq2 = 1.0/sqrt(2)
    sq6 = 1.0/sqrt(6)
    
    chik = [chi(sq6 * 2.0 * x, sq6 * x + sq2 * y, sq6 * x - sq2 * y) for x in k, y in k]
    
    fig, ax, hm = heatmap(k,k,chik,axis = (;aspect = 1))
    ax.title = "Cubic lattice T-flow"
    Colorbar(fig[1,2],hm)
    display("image/png", fig)
end