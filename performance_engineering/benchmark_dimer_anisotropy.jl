#!/usr/bin/env julia
# Benchmarking script for dimer anisotropy example using BenchmarkTools.jl
# Based on pmfrg_xyz/example-dimer-anisotropy.sh

using BenchmarkTools

# Activate parent project to access dependencies
using Pkg
Pkg.activate(dirname(@__DIR__))

using SpinFRGLattices
include("../Tpmfrg_xyz.jl")
using .Tpmfrg_xyz

println("Setting up dimer anisotropy system...")
system = SpinFRGLattices.getPolymer(2)
par = Params(system)
isotropy = zeros(3, length(system.couplings))
isotropy .= [1.0, 0.5, 0.2]
isotropy_matrix = Matrix(transpose(isotropy))

println("Running initial execution (for compilation)...")
SolveFRG(par, isotropy_matrix)

println("\nBenchmarking SolveFRG...")
b = @benchmark SolveFRG($par, $isotropy_matrix)

println("\nBenchmark results:")
display(b)
println("\n")
println("Minimum time: ", minimum(b.times) / 1e9, " seconds")
println("Median time:  ", median(b.times) / 1e9, " seconds")
println("Mean time:    ", mean(b.times) / 1e9, " seconds")
println("Maximum time: ", maximum(b.times) / 1e9, " seconds")
