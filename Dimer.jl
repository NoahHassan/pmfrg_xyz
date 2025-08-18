using Pkg
Pkg.activate(".")

include("Tpmfrg_xyz.jl")
using .Tpmfrg_xyz

using JLD2
using RecursiveArrayTools
using SpinFRGLattices,OrdinaryDiffEq,DiffEqCallbacks,RecursiveArrayTools,StructArrays
using SpinFRGLattices.StaticArrays
using CairoMakie

System = getPolymer(2)
isotropy = zeros(System.Npairs, 3)
for n in System.Npairs
    isotropy[n, :] = [0.4, 0.0, 1.0] # [0.0, 1.0, 0.0]
end

let
    Par = Tpmfrg_xyz.Params(
        System,
        N = 8,
        accuracy = 1e-10,
        temp_max = 100.,
        temp_min = 1.0
    )

    sol, saved_values = Tpmfrg_xyz.SolveFRG(Par, isotropy, method = Euler());
    save_object("ChiDimer.jld2", [(saved_values.saveval[n], exp(saved_values.t[n])) for n in 1:length(saved_values.t)])
end

function CompareArrays(arr1, arr2)
    mismatches = [Tuple(I) for I in CartesianIndices(arr1) if arr1[I] != arr2[I]]
    for n in 1:length(mismatches)
        println(mismatches[n])
    end
end

g1_before = load_object("tests_1/X_1.jld2")
g1_after  = load_object("tests_1/X_after_1.jld2")
g2_before = load_object("tests_2/X_1.jld2")
g2_after  = load_object("tests_2/X_after_1.jld2")

### Test iSigma
sig1_before = g1_before[1].iSigma
sig1_after  = g1_after[1].iSigma
sig2_before = g2_before[1].iSigma
sig2_after  = g2_after[1].iSigma

CompareArrays(sig1_after.x, sig2_before.x)
CompareArrays(sig1_after.y, sig2_before.y)
CompareArrays(sig1_after.z, sig2_before.z)

### Test Gamma
gam1_before = g1_before[1].Gamma
gam1_after  = g1_after[1].Gamma
gam2_before = g2_before[1].Gamma
gam2_after  = g2_after[1].Gamma

CompareArrays(gam1_after, gam2_before)

### Test X
X1_before = g1_before[2]
X1_after  = g1_after[2]
X2_before = g2_before[2]
X2_after  = g2_after[2]

CompareArrays(X1_after, X2_before)