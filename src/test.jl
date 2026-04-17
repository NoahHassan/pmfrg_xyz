using Pkg
Pkg.activate(".")

include("PMFRG_general.jl")
using .PMFRG_general

include("PMFRG_xyz.jl")
using .PMFRG_xyz

include("PMFRG_old.jl")
using .PMFRG_old

using SpinFRGLattices
using JLD2

System = getPolymer(2)

anisotropy = zeros(System.Npairs, 3, 3)
for Rij in 1:System.Npairs
    anisotropy[Rij, :, :] = [
        -0.2 0.0 0.0;
        0.0 3.0 0.0;
        0.0 0.0 1.0
    ]
end

Par = PMFRG_general.Params(
    System,
    N=4,
    temp_max=100.0,
    temp_min=1.0,
    accuracy=1e-6
)

PMFRG_general.TestFRG(Par, anisotropy)

let
    sol, saved_values = PMFRG_general.SolveFRG(Par, anisotropy)
    save_object("general_[-0.2 3.0 1.0]_polymer.jld2", [(saved_values.saveval[n], exp(saved_values.t[n])) for n in 1:length(saved_values.t)])
end

System = getPolymer(2)
anisotropy = zeros(System.Npairs, 3)
for Rij in 1:System.Npairs
    anisotropy[Rij, :] = [-0.2, 3.0, 1.0]
end

Par = PMFRG_xyz.Params(
    System,
    N=4,
    temp_max=100.0,
    temp_min=1.0,
    accuracy=1e-6
)

Par_old = PMFRG_old.Params(
    System,
    N=4,
    temp_max=100.0,
    temp_min=1.0,
    accuracy=1e-6
)

PMFRG_xyz.TestFRG(Par, anisotropy)
PMFRG_old.TestFRG(Par_old, anisotropy)

let
    sol, saved_values = PMFRG_xyz.SolveFRG(Par, anisotropy)
    save_object("xyz_[-0.2 3.0 1.0]_polymer.jld2", [(saved_values.saveval[n], exp(saved_values.t[n])) for n in 1:length(saved_values.t)])
end

d_gen = load_object("general_[-0.2 3.0 1.0]_polymer.jld2")
d_xyz = load_object("xyz_[-0.2 3.0 1.0]_polymer.jld2")

chi_gen = getindex.(d_gen, 1)
chi_xyz = getindex.(d_xyz, 1)
Tvals = getindex.(d_gen, 2)

let
    n=600
    println(chi_gen[n].Chi_z .- chi_xyz[n].Chi_z)
end
#0.0024367562701372576

a = -5.937781120053247e-6
b = -1.187556224010651e-5
2 * a

# only adding GGGG into chi_z for the values at which Γ is not zero gives
a = 2.37511244802127e-5 # general
b = -5.937781120053247e-6 # xyz
# the general value is 4 times as large as the xyz value
# since i output *all* Γ that are non-zero and they are the same for xyz and general
# the error must be with GGGG
# since otherwise the error is so consistent, the flow equations are probably even correct...

# TestFRG: 
# chi_gen
# [0.0024367562701372576, 1.187556224010646e-6]
# [0.0024367562701372576, -1.781334336015969e-5]
# [0.0024367562701372576, -5.937781120053247e-6]
#
# chi_xyz
# [0.0024367562701372576, 1.187556224010646e-6]
# [0.0024367562701372576, -1.781334336015969e-5]
# [0.0024367562701372576, -5.937781120053247e-6]
#
# Hence first getChi_z must be correct.

# Monitor now iSigma and Gamma throughout the entire flow

using DiffEqCallbacks

struct Observables_Full{T}
    Chi_x::Vector{T}
    Chi_y::Vector{T}
    Chi_z::Vector{T}
    iSigma::Array{T, 3} # 9 Sigma-flavors
    Gamma::Array{T, 5} # 81 Gamma-flavors
end

function save_all(State, t, Par)
    chix = PMFRG_general.getChi_x(State, exp(t), Par)
    chiy = PMFRG_general.getChi_y(State, exp(t), Par)
    chiz = PMFRG_general.getChi_z(State, exp(t), Par)
    println("t: $(round(exp(t);digits=3))")
    return Observables_Full(
        copy(chix), copy(chiy), copy(chiz),
        copy(State.x[2]), copy(State.x[3])
    )
end

System = getPolymer(2)

anisotropy = zeros(System.Npairs, 3, 3)
for Rij in 1:System.Npairs
    anisotropy[Rij, :, :] = [
        -0.2 0.0 0.0;
        0.0 3.0 0.0;
        0.0 0.0 1.0
    ]
end

Par = PMFRG_general.Params(
    System,
    N=4,
    temp_max=100.0,
    temp_min=1.0,
    accuracy=1e-6
)

sol, saved_values =
    PMFRG_general.SolveFRG(
        Par, anisotropy,
        SavedValues(PMFRG_general._getFloatType(Par), Observables_Full{PMFRG_general._getFloatType(Par)}),
        (State, t, integrator) -> save_all(State, t, Par)
    )

save_object("general_[-0.2 3.0 1.0]_polymer.jld2", [(saved_values.saveval[n], exp(saved_values.t[n])) for n in 1:length(saved_values.t)])

struct Observables_Full_XYZ{T}
    Chi_x::Vector{T}
    Chi_y::Vector{T}
    Chi_z::Vector{T}
    iSigma_x::Array{T, 2}
    iSigma_y::Array{T, 2}
    iSigma_z::Array{T, 2}
    Gamma::Array{T, 5}
end

function save_all_xyz(State, t, Par)
    chix = PMFRG_xyz.getChi_x(State, exp(t), Par)
    chiy = PMFRG_xyz.getChi_y(State, exp(t), Par)
    chiz = PMFRG_xyz.getChi_z(State, exp(t), Par)
    println("t: $(round(exp(t); digits=3))")
    return Observables_Full_XYZ(
        copy(chix), copy(chiy), copy(chiz),
        copy(State.x[2]), copy(State.x[3]), copy(State.x[4]),
        copy(State.x[5])
    )
end

System = getPolymer(2)

anisotropy = zeros(System.Npairs, 3)
for Rij in 1:System.Npairs
    anisotropy[Rij, :] = [-0.2, 3.0, 1.0]
end

Par = PMFRG_xyz.Params(
    System,
    N=4,
    temp_max=100.0,
    temp_min=1.0,
    accuracy=1e-6
)

let
    sol, saved_values =
        PMFRG_xyz.SolveFRG(
            Par,
            anisotropy,
            SavedValues(PMFRG_xyz._getFloatType(Par), Observables_Full_XYZ{PMFRG_xyz._getFloatType(Par)}),
            (State, t, integrator) -> save_all_xyz(State, t, Par)
        )
    save_object("xyz_[-0.2 3.0 1.0]_polymer.jld2", [(saved_values.saveval[n], exp(saved_values.t[n])) for n in 1:length(saved_values.t)])
end

data_gen = load_object("general_[-0.2 3.0 1.0]_polymer.jld2")
data_xyz = load_object("xyz_[-0.2 3.0 1.0]_polymer.jld2")

sig_gen = getfield.(getindex.(data_gen, 1), :iSigma)
sigx_xyz = getfield.(getindex.(data_xyz, 1), :iSigma_x)

gam_gen = getfield.(getindex.(data_gen, 1), :Gamma)
gam_xyz = getfield.(getindex.(data_xyz, 1), :Gamma)

# Before switching sign in addTo1PartBubble (n=1):
# [5.096792215915629e-9 1.7645491252554948e-9 1.1042344391208137e-9 7.171418595963982e-10]
# [0.0, -1.3111058764277844e-8, 0.0, -1.2991740763648296e-8]
# [0.0, -1.2210483735231037e-9, 0.0, -1.7117963713175754e-9]
#
# After switching sign in addTo1PartBubble (n=2):
# [5.096792215915629e-9 1.7645491252554948e-9 1.1042344391208137e-9 7.171418595963982e-10]
# [0.0, -1.3111058764277844e-8, 0.0, -1.2991740763648296e-8]
# [0.0, -1.2210483735231037e-9, 0.0, -1.7117963713175754e-9]

let
    n = 400
    println(sig_gen[n][1, :, :] .- sigx_xyz[n])
    println(gam_gen[n][1, 2, 1, 1, :] .- gam_xyz[n][1, 2, 1, 1, :]) # xxxx
    println(gam_gen[n][11, 2, 1, 1, :] .- gam_xyz[n][10, 2, 1, 1, :]) # xyxy
end

x_gen = load_object("X_gen_test.jld2")
x_xyz = load_object("X_xyz_test.jld2")

y_gen = load_object("Y_gen.jld2")
y_xyz = load_object("Y_xyz.jld2")

for n in eachindex(y_xyz)
    if(y_xyz[n] != 0.0)
        println("$n: $(y_xyz[n])")
    end
end

function IndexFromFlavors(d1::Int, d2::Int, d3::Int, d4::Int)
    return d1 * 27 + d2 * 9 + d3 * 3 + d4 + 1
end

y_gen[IndexFromFlavors(2,1,2,1)]

33-21

# The values of Y_xyz and Y_gen are identical in the first step

g_gen = load_object("Gamma_gen.jld2")[:, 2, 1, 1, :]
g_xyz = load_object("Gamma_xyz.jld2")[:, 2, 1, 1, :]

println(g_xyz[5, :])
println(g_gen[IndexFromFlavors(0,0,2,2), :])

# Dgammas are equal up to the second step

dg_gen = load_object("Dgamma_gen.jld2")
dg_xyz = load_object("Dgamma_xyz.jld2")

dg_gen[9, :, :]
dg_xyz.z