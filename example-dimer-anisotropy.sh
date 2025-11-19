#!/bin/bash

time julia --project=. -e '
using SpinFRGLattices
include("Tpmfrg_xyz.jl")
using .Tpmfrg_xyz
system = SpinFRGLattices.getPolymer(2)
par = Params(system)
isotropy = zeros(3,length(system.couplings))
isotropy .= [1.0,0.5,0.2]
SolveFRG(par,Matrix(transpose(isotropy)))
'
