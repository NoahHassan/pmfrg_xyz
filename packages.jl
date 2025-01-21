### pinned SpinFRGLattices at v0.5.2
import Pkg;

### so that the debugger doesnt always execute this code
if !@isdefined(pkgs_loaded)
    println("Loading packages...")
    Pkg.add("SpinFRGLattices")

    Pkg.add("OrdinaryDiffEq")
    Pkg.add("RecursiveArrayTools")
    Pkg.add("CairoMakie")
    Pkg.add("DiffEqCallbacks")
    Pkg.add("StructArrays")
    Pkg.add("JLD2")

    pkgs_loaded = true
end