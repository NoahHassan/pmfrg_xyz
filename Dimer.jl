using Pkg
Pkg.activate(".")

include("Tpmfrg_xyz.jl")
using .Tpmfrg_xyz

include("pmfrg_xyz.jl")
using .pmfrg_xyz

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
        temp_min = 0.1
    )

    sol, saved_values = Tpmfrg_xyz.SolveFRG(Par, isotropy, method = DP5());
    save_object("ChiDimer_TFlow_Low.jld2", [(saved_values.saveval[n], exp(saved_values.t[n])) for n in 1:length(saved_values.t)])
end

# Test with Λ-Flow
let
    Trange = LinRange(1.0, 0.1, 100)

    chis = []
    for T_ in Trange
        println("T = $T_")
        Par = pmfrg_xyz.Params(
            System,
            T = T_,
            N = 8,
            accuracy = 1e-10,
            lambda_max = 100.,
            lambda_min = 0.1
        )

        sol, saved_values = pmfrg_xyz.SolveFRG(Par, isotropy, method = DP5());
        append!(chis, (saved_values.saveval[end], T_))
    end
    save_object("ChiDimer_LFlow_Low.jld2", chis)
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