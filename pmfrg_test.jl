using Pkg
Pkg.activate(".")
include("Tpmfrg_xyz.jl")
import .Tpmfrg_xyz

using JLD2
using SpinFRGLattices, OrdinaryDiffEq, DiffEqCallbacks, RecursiveArrayTools, StructArrays
using SpinFRGLattices.StaticArrays
using SpinFRGLattices.SquareLattice

J1 = 1.0
J2 = 0.5

System = getSquareLattice(6, [J1, J2])
isotropy = zeros(System.Npairs, 3)

for n in 1:System.Npairs
    isotropy[n, :] = [1.0, 0.5, 0.2]
end

Par = Tpmfrg_xyz.Params(
    System,
    N = 8,
    accuracy = 1e-10,
    temp_max = 10.0,
    temp_min = 1.0
)

Tpmfrg_xyz.SolveFRG(Par, isotropy, method=DP5())